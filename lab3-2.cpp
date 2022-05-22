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
#include <list>
#include <random>
static std::default_random_engine engine(10) ; // random s e e d = 10
static std::uniform_real_distribution<double> uniform ( 0 , 1 );

std::vector<int> ob_list;

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
Vector operator*(const Vector &a, const Vector &b)
{
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
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

class Geometry
{
public:
	Geometry(){};
	virtual bool intersection(const Ray& ray, Vector &P, Vector &N, double &t,
                              const Vector& A, const Vector& B, const Vector& C, 
                              double &alpha, double &beta, double &gamma, int &id_mesh) = 0;
	Vector albedo;
	bool isMirror,isTrans;
};


class Sphere: public Geometry {
public:
	Vector center;
	double r;

	Sphere(const Vector& center, double r, const Vector& color,bool ism = false, bool ist = false):Geometry(){
        this->center=center;
        this->r=r;
        this->albedo=color;
        this->isMirror=ism;
        this->isTrans=ist;
    }

	bool intersection(const Ray& ray, Vector &P, Vector &N, double &t,
                              const Vector& A, const Vector& B, const Vector& C, 
                              double &alpha, double &beta, double &gamma, int &id_mesh){
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

Vector random_cos(const Vector &N){
    double r1,r2;
    r1 = uniform(engine);
    r2 = uniform(engine);

    double x,y,z;
    x = cos(2*M_PI* r1) * sqrt(1 - r2);
    y = sin(2*M_PI* r1) * sqrt(1 - r2);
    z = sqrt(r2);

    Vector T1,T2;
    double min_num = N[0];
    if (N[1]<min_num){
        min_num = N[1];
    }
    if (N[2]<min_num){
        min_num = N[2];
    }

    if (min_num == N[0]){
        T1 = Vector(0, N[2], -N[1]);
    }
    if (min_num == N[1]){
        T1 = Vector(N[1], 0, -N[0]);
    }
    if (min_num == N[2]){
        T1 = Vector(N[1], -N[0], 0);
    }
    T2 = cross(N, T1);

    T1.normalize();
    T2.normalize();
    return x*T1 + y*T2 + z*N;
}

class BoundingBox 
{
public:
    Vector min,max;
    BoundingBox(Vector min = Vector(0, 0, 0), Vector max = Vector(0, 0, 0)) : min(min), max(max) {};
    bool intersect(const Ray &ray){
        double t0[3],t1[3];
        for (int i = 0; i<3; i++){
            double x,y;
            x = std::min( (min[i] - ray.O[i]) / ray.dirrection[i], (max[i] - ray.O[i]) / ray.dirrection[i] );
            y = std::max( (min[i] - ray.O[i]) / ray.dirrection[i], (max[i] - ray.O[i]) / ray.dirrection[i] );
            t0[i] = x;
            t1[i] = y;
        }

        double maxt0, mint1;
        maxt0 = std::max(t0[0],std::max(t0[1],t0[2]));
        mint1 = std::min(t1[0],std::min(t1[1],t1[2]));

        if (mint1 > maxt0 and maxt0 > 0){
            return true;
        }

        return false;
        
    }
};

void boxMuller ( double stdev , double &x , double &y ) {
	double r1 = uniform ( engine ) ;
	double r2 = uniform ( engine ) ;
	x = sqrt(-2 * log ( r1 ) ) *cos ( 2 * M_PI*r2 ) *stdev ;
	y = sqrt(-2 * log ( r1 ) ) *sin ( 2 * M_PI*r2 ) *stdev ;
}

bool Triangle_Intersection (const Ray& ray, const Vector& A, const Vector& B, const Vector& C, double &alpha, double &beta, double &gamma, double &t, Vector &N){
	Vector e1,e2;
	e1 = B-A;
	e2 = C-A;
	N = cross(e1,e2);
	double b1,y1,dv;
	b1 = dot(e2,cross((A - ray.O),ray.dirrection));
	y1 = dot(e1,cross((A - ray.O),ray.dirrection));
	dv = dot(ray.dirrection,N);

	beta = b1/dv;
	gamma = -y1/dv;
	alpha = 1-beta-gamma;
	t = dot(A-ray.O,N)/dot(ray.dirrection,N);

	if (alpha<0 or alpha>1 or beta<0 or beta>1 or gamma<0 or gamma>1 or t<0){
		return false;
	}
	return true;
}

class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};

class BVH{
public:
    int head,tail;
    BoundingBox box;
    BVH *left, *right;
};

class TriangleMesh: public Geometry {
public:
    ~TriangleMesh() {}
	TriangleMesh(const Vector color) {
        this->albedo = color;
		this->isMirror = false;
		this->isTrans = false;
    };

	bool intersection(const Ray& ray, Vector &P, Vector &N, double &t,
                      const Vector& A, const Vector& B, const Vector& C, 
                      double &alpha, double &beta, double &gamma,int &id_mesh){
        //std::cout << "GG " << std::endl;
        if (! root.box.intersect(ray)){
            return false;
        }

        std::list<BVH*> nodes_to_visit;
        nodes_to_visit.push_front(&root);
        while (! nodes_to_visit.empty()){
            BVH* curNode = nodes_to_visit.back();
            nodes_to_visit.pop_back();

            if (curNode->left != NULL and curNode->right != NULL){
                if (curNode->left->box.intersect(ray)){
                    nodes_to_visit.push_back(curNode->left);
                }

                if (curNode->right->box.intersect(ray)){
                    nodes_to_visit.push_back(curNode->right);
                }
            }
            if ((curNode->left == NULL) || (curNode->right == NULL)){
                t = 3E10;
                bool f_inter = false;
                int id = -1;
                //std::cout << "gg " << std::endl;
                //std::cout << "head tail " << curNode->head << "| " << curNode->tail << std::endl;
                for (int i = curNode->head; i < curNode->tail; i++) {
                    //std::cout << "i " << i << std::endl;
                    Vector localN;
                    double localt, localalpha, localbeta, localgamma;
                    if ( Triangle_Intersection(ray, vertices[indices[i].vtxi], vertices[indices[i].vtxj], vertices[indices[i].vtxk], localalpha, localbeta, localgamma, localt, localN))
                    {
                        if (localt < t)
                        {
                            t = localt;
                            P = ray.O + ray.dirrection * t;
                            N = localN;
                            id = i;
                            f_inter = true;
                        }
                    }
                }
                id_mesh = id;
                return f_inter;
            }
        }
	}

    void Compute_min_max(Vector &min, Vector &max,int l, int r){
        min = vertices[indices[l].vtxi];
        max = vertices[indices[l].vtxi];
        for (int i = l; i < r; i++){
            min[0] = std::min(min[0],vertices[indices[i].vtxi][0]);
            min[1] = std::min(min[1],vertices[indices[i].vtxi][1]);
            min[2] = std::min(min[2],vertices[indices[i].vtxi][2]);
            max[0] = std::max(max[0],vertices[indices[i].vtxi][0]);
            max[1] = std::max(max[1],vertices[indices[i].vtxi][1]);
            max[2] = std::max(max[2],vertices[indices[i].vtxi][2]);

            min[0] = std::min(min[0],vertices[indices[i].vtxj][0]);
            min[1] = std::min(min[1],vertices[indices[i].vtxj][1]);
            min[2] = std::min(min[2],vertices[indices[i].vtxj][2]);
            max[0] = std::max(max[0],vertices[indices[i].vtxj][0]);
            max[1] = std::max(max[1],vertices[indices[i].vtxj][1]);
            max[2] = std::max(max[2],vertices[indices[i].vtxj][2]);

            min[0] = std::min(min[0],vertices[indices[i].vtxk][0]);
            min[1] = std::min(min[1],vertices[indices[i].vtxk][1]);
            min[2] = std::min(min[2],vertices[indices[i].vtxk][2]);
            max[0] = std::max(max[0],vertices[indices[i].vtxk][0]);
            max[1] = std::max(max[1],vertices[indices[i].vtxk][1]);
            max[2] = std::max(max[2],vertices[indices[i].vtxk][2]);
        }
    }

    void recursive_call(BVH *H, int l, int r){
        Vector min,max;
        Compute_min_max(min,max,l,r);
        H->box = BoundingBox(min,max);
        H->head = l;
        H->tail = r;
        H->left = NULL;
        H->right = NULL;
        Vector diag = max - min;
        int diag_max = 0;
        if (diag[1] > diag[0]){
            diag_max = 1;
        }
        if (diag[2] > diag[diag_max]){
            diag_max = 2;
        }

        int pivot_index = l;

        Vector middle = diag*0.5 + min;
        double middle_axis = middle[diag_max];
        for (int i = l; i<r; i++){
            Vector barycenter = (vertices[indices[i].vtxi] + vertices[indices[i].vtxj] + vertices[indices[i].vtxk]) / 3;
            if (barycenter[diag_max] < middle_axis){
                std::swap ( indices [ i ] , indices [ pivot_index ] ) ;
                pivot_index++;
            }
        }

        if (pivot_index-l<=5 or r-pivot_index<=5 or r-l<=10){
            //std::cout << "pivot_index " << pivot_index << std::endl;
            return;
        }
        H->left = new BVH;
        H->right = new BVH;
        recursive_call ( H->left , l , pivot_index ) ;
        recursive_call ( H->right , pivot_index, r ) ;
   }

    void factor_move(double s, const Vector& t) {
        for (int i = 0; i < vertices.size(); i++) {
            vertices[i] = vertices[i] * s + t;
        }
    }
	
	void readOBJ(const char* obj) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;

				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
					col[0] = std::min(1., std::max(0., col[0]));
					col[1] = std::min(1., std::max(0., col[1]));
					col[2] = std::min(1., std::max(0., col[2]));

					vertices.push_back(vec);
					vertexcolors.push_back(col);

				} else {
					sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
					indices.push_back(t);
				} else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					} else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						} else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					} else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						} else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;								
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							} else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								} else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}

		}
		fclose(f);

	}

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
    BVH root;
	
};


class Scene {
public:
	Scene() {};
	void addObject(Geometry &s){
		objects.push_back(&s);
	}
	
	std::vector<Geometry*> objects;

	bool intersection(const Ray& ray, Vector &P, Vector &N, double &t, int &ID, int &id_mesh){
		t = 3E10;
		bool f_inter = false;
		for (int i =0; i<objects.size(); i++){
			Vector localP,localN;
			double localt;

            Vector A,B,C;
            double alpha,beta,gamma;
			id_mesh = -1;
			int localid_mesh;

			bool inter = objects[i]->intersection(ray, localP, localN, localt,A, B, C, alpha, beta, gamma,localid_mesh);
			if (inter and localt >0){
				f_inter = true;
				if (localt < t){
					ID = i;
					t = localt;
					P = localP;
					N = localN;
					id_mesh = localid_mesh;
				}
			}
		}

		return f_inter;
	}

	Vector getColor(const Ray& ray, int bounce){
		Vector color(0,0,0);
		Vector P,N;
		Vector L (-10,20,40);
		int ID;
		double t;
		double I = 2E10;
		int id_mesh;
		
		//std::cout << "intersection " << intersection(ray, P, N, t, ID) << std::endl;

		if (bounce <= 0){
			return color;
		}

		if (intersection(ray, P, N, t, ID, id_mesh)){
			if (objects[ID]->isMirror){
				Vector R = ray.dirrection - 2*dot(ray.dirrection,N) *N;
				Ray reflect_ray(P+0.001*N, R);
				return getColor(reflect_ray,bounce-1);
			}

			if (objects[ID]->isTrans){
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
			int IDlight, id_meshlight;
			double tlight;
			Ray lightRay(P+0.001*N, lightVec);
			double Shadow = 1.;
			if (intersection(lightRay,Plight,Nlight,tlight,IDlight,id_meshlight)){
				if (tlight*tlight < distlight){
					Shadow = 0;
				}
			}
			Vector L0 = Shadow * I/((L-P).norm2()*4*M_PI) * objects[ID]->albedo / M_PI * std::max(0.,dot(lightVec,N));  
            Ray ray_random(P, random_cos(N)); 
			color = L0 + objects[ID]->albedo * getColor(ray_random, bounce - 1);
		}
		return color;
	}
	
};


int main() {
    auto start = std::chrono::steady_clock::now();
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
	scene.addObject(Sfloor);
	scene.addObject(Sceilling);
	scene.addObject(Sleft);
	scene.addObject(Sright);
	scene.addObject(Sback);
	scene.addObject(Sfront);
    
	//std::cout << "addmesh" <<  std::endl;

    TriangleMesh mesh(Vector(0.3,0.3,0.3));
    mesh.readOBJ("cat.obj");
	mesh.factor_move(0.6, Vector(0, -10, 0));	
    mesh.recursive_call(&mesh.root, 0 ,mesh.indices.size());

    std::cout << mesh.indices.size() << std::endl;
    std::cout << "mesh vertices" << mesh.vertices.size() << std::endl;
    //std::cout << "mesh max min" << mesh.root.box.max[0] << mesh.root.box.max[1] << mesh.root.box.max[2] << "| " << mesh.root.box.min[0] << mesh.root.box.min[1] << mesh.root.box.min[2] << std::endl;

    scene.addObject(mesh);

	ob_list.push_back(1);
	ob_list.push_back(1);
	ob_list.push_back(1);
	ob_list.push_back(1);
	ob_list.push_back(1);
	ob_list.push_back(1);
	ob_list.push_back(2);

	//std::cout << "size: " << scene.objects.size() << '|' << ob_list.size() << std::endl;

	std::vector<unsigned char> image(W * H * 3, 0);
    int nb_paths = 20;
#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
            Vector color(0.,0.,0.);
            for (int k = 0; k < nb_paths; k++){

                double x, y;
                boxMuller(1, x, y);
                x = i;
                y = j;
                Vector u(y - W/2 + 0.5, H/2 - x - 0.5, -W / (2 * tan( fov / 2)));
                u.normalize();
                Ray ray(Q,u);
                color = color + scene.getColor(ray,5);   
            }
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0]/nb_paths, 1. / 2.2));
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1]/nb_paths, 1. / 2.2));
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2]/nb_paths, 1. / 2.2));
        
        }
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

    auto finish = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
     std::cout << "Time is " << elapsed << " microseconds" << std::endl;

	return 0;
}