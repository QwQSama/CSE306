#include <vector>
#include <iostream>
#include <math.h>

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

// if the Polygon class name conflicts with a class in wingdi.h on Windows, use a namespace or change the name
class Polygon {  
public:
	std::vector<Vector> vertices;
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

void lineIntersectSeg(Vector A, Vector B, Vector l1, Vector l2, Vector &point, bool &f) {
	Vector N = Vector(l2[1] - l1[1], -l2[0] + l1[0]);
	Vector u = l1;
	double t = dot(u - A, N) / dot(B - A, N);
	point = A + t * (B - A);
	f = true;
	if (((dot(u - B, N) <= 0) and (dot(u - A, N) > 0)) or ((dot(u - B, N) > 0) and (dot(u - A, N) <= 0))){
		f = false;
	}
	
}

Polygon clip_cell_edge(Polygon cell, Vector x, Vector l1, Vector l2){
    Polygon res;
    int n = cell.vertices.size();
	std::vector<Vector> position;
	std::vector<bool> inside;

	for (int i=0; i<n; i++){
		bool f;
		Vector p;
		lineIntersectSeg(cell.vertices[i],x, l1,l2,p,f);
		inside.push_back(f);

		int next = i+1;
		if (i == n-1){
			next = 0;
		}
		lineIntersectSeg(cell.vertices[i],cell.vertices[next], l1,l2,p,f);
		position.push_back(p);
    }

	if (inside[0]){
		res.vertices.push_back(cell.vertices[0]);
	}

	for (int i=1; i<n-1; i++){
		if (inside[i]){
			if (!inside[i-1]){
				res.vertices.push_back(position[i-1]);
			}
			res.vertices.push_back(cell.vertices[i]);
		}
		else{
			if (inside[i-1]){
				Vector p = position[i-1];
				res.vertices.push_back(p);
			}
		}
	}

	int i = n-1;
	if (inside[i]){
		if (!inside[i-1]){
			res.vertices.push_back(position[i-1]);
		}
		res.vertices.push_back(cell.vertices[i]);
		if (!inside[0]){
			res.vertices.push_back(position[i]);
		}
	}
	else{
		if (inside[i-1]){
			Vector p = position[i-1];
			res.vertices.push_back(p);
		}
		if (inside[0]){
			res.vertices.push_back(position[i]);
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

Polygon add_cell(Vector* points, int n, int js){
    Polygon res;
	res.vertices.push_back(Vector(0,0,0));
	res.vertices.push_back(Vector(0,1,0));
	res.vertices.push_back(Vector(1,1,0));
	res.vertices.push_back(Vector(1,0,0));
	
    for (int i=0; i<n; i++){
        if (i==js){
            continue;
        }
		
        Vector line_begin = (points[i] + points[js])/2;
        Vector line_end = xuanzhuan_pi_2(points[i],line_begin);

		/*
		if (js==0 and i==1){
			std::cout << "line" << line_begin.data[0] << line_begin.data[1] << line_end.data[0] << line_end.data[1] << std::endl;
			std::cout << "point" << points[i].data[0] << points[i].data[1] << points[js].data[0] << points[js].data[1] << std::endl;

			res = clip_cell_edge(res,points[js],line_begin,line_end);
		}
		/*
        std::cout << "dot" << dot(points[i]-points[js],line_begin-line_end) << std::endl;
		Vector p;
		lineIntersectSeg(Vector(0,0,0),Vector(0,1,0),Vector(1,0,0),Vector(0,0,0),p);
        std::cout << "intersec" << p.data[0] << p.data[1] << std::endl;
		*/

        res = clip_cell_edge(res,points[js],line_begin,line_end);
   }

	
	return res;
}

std::vector<Polygon> voronoi_diagram(int num_points, Vector* points){
	std::vector<Polygon> res;
	for (int i=0; i<num_points; i++){
		Polygon pyg = add_cell(&points[0],num_points,i);
		res.push_back(pyg);
	}
	return res;
}

int main(){
	
	int num_points = 400;
	std::vector<Vector> points(num_points);
	for (int i = 0; i < num_points; i++) 
	{
		for (int j = 0; j < 2; j++)
		{
			points[i][j] = rand() / (double)RAND_MAX;
		}
		points[i][2] = 0;
	}

	std::vector<Polygon> diagram = voronoi_diagram(num_points, &points[0]);
	save_svg(diagram, "diagram.svg");

    return 0;
}