#include <vector>
#include <iostream>
#include <math.h>
#include "lbfgs.h"

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
double cross(const Vector& a, const Vector& b) {
	//return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
	return a[0] * b[1] - a[1] * b[0];
}

// if the Polygon class name conflicts with a class in wingdi.h on Windows, use a namespace or change the name
class Polygon {  
public:
	std::vector<Vector> vertices;

	double area() {
		double a = 0;
		if (vertices.size() < 2) return 0;
		for (int i = 0; i < vertices.size(); i++) {
			int j = (i + 1) % vertices.size();
			a += (vertices[i][0] * vertices[j][1] - vertices[j][0] * vertices[i][1]);
		}
		return std::abs(0.5 * a);
	}

	double integral_squared_dist(Vector& P) {
		double a = 0;
		if (vertices.size() < 3) return 0;
		for (int i = 1; i < vertices.size() - 1; i++) {
			double areaT = 0.5 * cross(vertices[i] - vertices[0], vertices[i + 1] - vertices[0]);
			areaT = std::abs(areaT);
			Vector C[3] = {vertices[0], vertices[i], vertices[i + 1]};
			for (int l = 0; l < 3; l++) {
				for (int m = l; m < 3; m++) {
					a += areaT / 6 * dot(C[l] - P, C[m] - P);
				}

			}
		}
		
		return a;
	}
};	

// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
	void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
		FILE* f = fopen(filename.c_str(), "w+"); 
		fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        for (int i=0; i<polygons.size(); i++) {
		    fprintf(f, "<g>\n");
		    fprintf(f, "<polygon points = \""); 
		    for (int j = 0; j < polygons[i].vertices.size(); j++) {
			    fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
		    }
		    fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
		    fprintf(f, "</g>\n");
        }
		fprintf(f, "</svg>\n");
		fclose(f);
	}


// Adds one frame of an animated svg file. frameid is the frame number (between 0 and nbframes-1).
// polygons is a list of polygons, describing the current frame.
// The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
	void save_svg_animated(const std::vector<Polygon> &polygons, std::string filename, int frameid, int nbframes) {
		FILE* f;
		if (frameid == 0) {
			f = fopen(filename.c_str(), "w+");
			fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
			fprintf(f, "<g>\n");
		} else {
			f = fopen(filename.c_str(), "a+");
		}
		fprintf(f, "<g>\n");
		for (int i = 0; i < polygons.size(); i++) {
			fprintf(f, "<polygon points = \""); 
			for (int j = 0; j < polygons[i].vertices.size(); j++) {
				fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000-polygons[i].vertices[j][1] * 1000));
			}
			fprintf(f, "\"\nfill = \"none\" stroke = \"black\"/>\n");
		}
		fprintf(f, "<animate\n");
		fprintf(f, "	id = \"frame%u\"\n", frameid);
		fprintf(f, "	attributeName = \"display\"\n");
		fprintf(f, "	values = \"");
		for (int j = 0; j < nbframes; j++) {
			if (frameid == j) {
				fprintf(f, "inline");
			} else {
				fprintf(f, "none");
			}
			fprintf(f, ";");
		}
		fprintf(f, "none\"\n	keyTimes = \"");
		for (int j = 0; j < nbframes; j++) {
			fprintf(f, "%2.3f", j / (double)(nbframes));
			fprintf(f, ";");
		}
		fprintf(f, "1\"\n	dur = \"5s\"\n");
		fprintf(f, "	begin = \"0s\"\n");
		fprintf(f, "	repeatCount = \"indefinite\"/>\n");
		fprintf(f, "</g>\n");
		if (frameid == nbframes - 1) {
			fprintf(f, "</g>\n");
			fprintf(f, "</svg>\n");
		}
		fclose(f);
	}

void lineIntersectSeg(Vector A, Vector B, Vector l1, Vector l2, Vector &point, int &f,double w1, double w2) {
	Vector N = l2 - l1;
	Vector u = (l1 + l2)/2 + (w1-w2) / (2*(l2-l1).norm2()) * (l2-l1);
	double t = dot(u - A, N) / dot(B - A, N);
	point = A + t * (B - A);

	if ( ( (B-l1).norm2() - w1) <= ( (B-l2).norm2() - w2) ){
		if (((A - l1).norm2() - w1) > ((A - l2).norm2() - w2)){
			f=1;
		}
		f=0;
	}
	else{
		if (((A - l1).norm2() - w1) <= ((A - l2).norm2() - w2)){
			f=-1;
		}
		f = -2;
	}
	
	
}

Polygon clip_cell_edge(Polygon cell, Vector* points, int js1, int js2, double w1 = 0.0, double w2 = 0.0){
	
    Polygon res;
    int n = cell.vertices.size();
	std::vector<Vector> position;
	std::vector<int> difside;

	for (int i=0; i<n; i++){
		int f;
		Vector p;

		int next = i+1;
		if (i == n-1){
			next = 0;
		}
		Vector A,B,l1,l2;
		A = cell.vertices[i];
		B = cell.vertices[next];
		l1 = points[js1];
		l2 = points[js2];
		
		Vector N = l2 - l1;
		Vector u = (l1 + l2)/2 + (w1-w2) / (2*(l2-l1).norm2()) * (l2-l1);
		double t = dot(u - A, N) / dot(B - A, N);
		Vector	point = A + t * (B - A);

		if ( ( (B-l1).norm2() - w1) <= ( (B-l2).norm2() - w2) ){
			if (((A - l1).norm2() - w1) > ((A - l2).norm2() - w2)){
				res.vertices.push_back(point);
			}
			res.vertices.push_back(B);
		}
		else{
			if (((A - l1).norm2() - w1) <= ((A - l2).norm2() - w2)){
				res.vertices.push_back(point);
			}
		}
	}

	return res;
}

Vector xuanzhuan_pi_2(Vector p, Vector p0){
	double x,y,dx,dy;
	x = p.data[0];
	y = p.data[1];
	dx = p0.data[0];
	dy = p0.data[1];
	double xx = -(y-dy) + dx;
	double yy = (x-dx) + dy;
	return Vector(xx,yy,0);
	
}

Polygon add_cell(Vector* points, int n, int js, const double* w){
    Polygon res;
	res.vertices.push_back(Vector(0,0,0));
	res.vertices.push_back(Vector(0,1,0));
	res.vertices.push_back(Vector(1,1,0));
	res.vertices.push_back(Vector(1,0,0));
	
    for (int i=0; i<n; i++){
        if (i==js){
            continue;
        }
		
		double w1, w2;
		w1 = w[js];
		w2 = w[i];
        //Vector line_begin = (points[i] + points[js])/2;
        //Vector line_end = xuanzhuan_pi_2(points[i],line_begin);

        res = clip_cell_edge(res,&points[0],js,i,w1,w2);
   }

	
	return res;
}

std::vector<Polygon> voronoi_diagram(int num_points, Vector* points, const double* w){
	std::vector<Polygon> res;
	for (int i=0; i<num_points; i++){
		Polygon pyg = add_cell(&points[0],num_points,i,&w[0]);
		res.push_back(pyg);
	}
	return res;
}

	
static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
    )
{
    int i;
    lbfgsfloatval_t fx = 0.0;

	Vector* points = static_cast<Vector*>(instance);
	std::vector<Polygon> diagram(n);
	diagram = voronoi_diagram(n, &points[0], &x[0]);

	double lambda = 1./n;
    for (i = 0; i < n; i += 1) {
		double area = diagram[i].area();
        g[i] = -(lambda - area);
		fx += -(diagram[i].integral_squared_dist(points[i]) - x[i] * area + lambda * x[i]);
	}

    return fx;
}

int main(){
	
	int num_points = 300;
	std::vector<Vector> points(num_points);
	std::vector<double> w(num_points);
	for (int i = 0; i < num_points; i++) 
	{
		for (int j = 0; j < 2; j++)
		{
			points[i][j] = rand() / (double)RAND_MAX;
		}
		points[i][2] = 0;

		w[i] = 0;
	}

	

	lbfgs_parameter_t param;
	lbfgs_parameter_init(&param);
    param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
	double fx;
	
	int ret = lbfgs(num_points, &w[0], &fx, evaluate, NULL, &points[0], &param);
	

	std::vector<Polygon> diagram = voronoi_diagram(num_points, &points[0],&w[0]);
	save_svg(diagram, "diagram.svg");

    return 0;
}