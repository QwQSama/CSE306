#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <math.h>

#define M_PI  3.14159265358979323846

#include <iostream>
#include <chrono>
#include <random>
static std::default_random_engine engine(10) ; // random s e e d = 10
static std::uniform_real_distribution<double> uniform ( 0 , 1 );

class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		data[0] = x;
		data[1] = y;
		data[2] = z;
	}
	double norm2() const {
		return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
	}
	double norm() const {
		return sqrt(norm2());
	}
	void normalize() {
		double n = norm();
		data[0] /= n;
		data[1] /= n;
		data[2] /= n;
	}
	double operator[](int i) const { return data[i]; };
	double& operator[](int i) { return data[i]; };
	double data[3];
};
static inline double sqr(double x) {
	return x*x;
}

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator-(const Vector& b) {
	return Vector( -b[0], -b[1], -b[2]);
}
Vector operator*(const double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

class Ray {
public:
	Vector O,dirrection;
	Ray(const Vector& pos, const Vector& dir):O(pos), dirrection(dir) {}
};

class Sphere {
public:
	Vector center;
	double r;
	Vector albedo;
	bool isMirror,isTrans;

	Sphere(const Vector& center, double r, const Vector& color,bool ism = false, bool ist = false): 
	center(center), r(r), albedo(color), isMirror(ism), isTrans(ist) {};

	bool intersection(const Ray& ray, Vector &P, Vector &N, double &t){
		double dotval = dot(ray.dirrection,ray.O - center);
		double delta = dotval*dotval - ((ray.O - center).norm2() - r*r);
		if (delta<0){
			return false;
		}

		double sqrtdelta = sqrt(delta);
		double t1 = -dotval - sqrtdelta;
		double t2 = -dotval + sqrtdelta;

		if (t2<0){
			return false;
		}

		if (t1>0){
			t = t1;
			P = ray.O + t1 * ray.dirrection;
		}
		else{
			t = t2;
			P = ray.O + t2 * ray.dirrection;
		}

		N = P - center;
		N.normalize();
		return true;
	}

};


class Scene {
public:
	Scene() {};
	void addObject(const Sphere& s){
		objects.push_back(s);
	}
	
	std::vector<Sphere> objects;

	bool intersection(const Ray& ray, Vector &P, Vector &N, double &t, int &ID){
		t = 2E10;
		bool f_inter = false;
		for (int i =0; i<objects.size(); i++){
			Vector localP,localN;
			double localt;

			bool inter = objects[i].intersection(ray,localP,localN,localt);
			if (inter){
				f_inter = true;
				if (localt < t){
					ID = i;
					t = localt;
					P = localP;
					N = localN;
				}
			}
		}

		return f_inter;
	}

	Vector getColor(const Ray& ray, int bounce){
		Vector color(0,0,0);
		Vector P,N;
		Vector L(-10,20,40);
		int ID;
		double t;
		double I = 1E10;

		if (bounce <= 0){
			return color;
		}

		if (intersection(ray, P, N, t, ID)){
			if (objects[ID].isMirror){
				Vector R = ray.dirrection - 2*dot(ray.dirrection,N) *N;
				Ray reflect_ray(P+0.001*N, R);
				return getColor(reflect_ray,bounce-1);
			}

			if (objects[ID].isTrans){
				Vector R = ray.dirrection - 2*dot(ray.dirrection,N) *N;
				Ray reflect_ray(P+0.001*N, R);

				double n1 =1;
				double n2 = 1.4;
				Vector N_Trans = N;
				if (dot(ray.dirrection, N) > 0){
					std::swap(n1,n2);
					N_Trans = -N_Trans;
					}

				double radic = 1-sqr(n1/n2) * (1-sqr(dot(ray.dirrection,N_Trans)));
				Vector Tt = n1/n2*(ray.dirrection - dot(ray.dirrection,N_Trans) *N_Trans);
				if (radic <0){
					return getColor(reflect_ray,bounce-1);
				}
				Vector Tn = -sqrt(radic)*N_Trans;
				Vector T = Tt + Tn;
				Ray refrect_ray (P+0.001*T, T);
				return getColor(refrect_ray,bounce-1);
			}

			Vector lightVec = (L-P);
			double distlight = lightVec.norm2();
			lightVec.normalize();

			Vector Plight,Nlight;
			int IDlight;
			double tlight;

			Ray lightRay(P+0.001*N, lightVec);
			double Shadow = 1.;
			if (intersection(lightRay,Plight,Nlight,tlight,IDlight)){
				if (tlight*tlight < distlight){
					Shadow = 0;
				}
			}
			color = Shadow * I/((L-P).norm2()*4*M_PI) * objects[ID].albedo / M_PI * std::max(0.,dot(lightVec,N));
		}
		return color;
	}
	
};

void boxMuller ( double stdev , double &x , double &y ) {
	double r1 = uniform ( engine ) ;
	double r2 = uniform ( engine ) ;
	x = sqrt(-2 * log ( r1 ) ) *cos ( 2 * M_PI*r2 ) *stdev ;
	y = sqrt(-2 * log ( r1 ) ) *sin ( 2 * M_PI*r2 ) *stdev ;
}



int main() {
	int W = 512;
	int H = 512;

	double fov = 60 * M_PI / 180;
	Vector Q(0,0,55);

	Sphere S1       (Vector(-20,0,0),    10, Vector(0.7,0.3,0.1));
	Sphere S2       (Vector(0,0,0),     10, Vector(0.7,0.3,0.1),true);
	Sphere S3       (Vector(20,0,0),     10, Vector(0.7,0.3,0.1),false,true);

	Sphere Sfloor   (Vector(0,-1000,0), 990, Vector(0.1,0.5,0.1));
	Sphere Sceilling(Vector(0,1000,0),  940, Vector(0.2,0.4,0.2));
	Sphere Sleft    (Vector(-1000,0,0), 940, Vector(0.1,0.3,0.8));
	Sphere Sright   (Vector(1000,0,0),  940, Vector(0.7,0.2,0.5));
	Sphere Sback    (Vector(0,0,1000),  940, Vector(0.1,0.4,0.2));
	Sphere Sfront   (Vector(0,0,-1000), 940, Vector(0.9,0.3,0.5));

	Scene scene;
	scene.addObject(S1);
	scene.addObject(S2);
	scene.addObject(S3);

	scene.addObject(Sfloor);
	scene.addObject(Sceilling);
	scene.addObject(Sleft);
	scene.addObject(Sright);
	scene.addObject(Sback);
	scene.addObject(Sfront);


	std::vector<unsigned char> image(W * H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			Vector u(j - W/2 + 0.5, H/2 - i - 0.5, -W / (2 * tan( fov / 2)));
			u.normalize();
			Ray ray(Q,u);
			Vector color = scene.getColor(ray,5);
			
			
			image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0], 0.45));
			image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1], 0.45));
			image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2], 0.45));
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}