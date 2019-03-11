#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//--------------Bounding Box  Class--------------------
class Bbox{
	public:
	double Xa[2];
	double Xb[2];
	void set_Xa(double x, double y);
	void set_Xb(double x, double y);
	void draw();
	void draw(char fn[128],char mode[3]);
	private:
};
Bbox bbox_union(Bbox b1, Bbox b2);
Bbox bbox_cross(Bbox b1, Bbox b2);
//------------------------Circ Class--------------------
class Circ{
	public:
		double xc[2];	// center 
		double radi;	// radius
		void draw(char fn[128],int npnt, char mode[3]);
		Bbox bbox;	// bounding box
		void set_bbox();	// set bounding box
		bool isin(double* x);
	private:
};
//------------------- Rectangle Class ------------------
class RECT{
	public:
	double x[4],y[4];
	double xmin,xmax,ymin, ymax;
	void setup(
		double xc[2], //	center
		double b,  //	width
		double h,  //	height
		double th  //	counterclock-wise angle in degree
	);
	void setup(
		double x1[2], //endpoint 1
		double x2[2], //endpoint 2
		double t,	// thickness
		int endcap	// endcap (0:none, 1:add)
	);
	void out(char *fname);
	int isin(double xf[2]);
	private:
};
class Ellip{
	public:
		double xc[2];	// center
		double radi[2]; // radii
		double phi;	// angle [rad]
		void draw(char fn[128], int np, char mode[3]);
		void draw(int np);
		void set_xc(double x,double y);
		void set_radi(double r1,double r2);
		void set_phi(double ang);
		void scale(double s);
		Bbox bbox;	// bounding box
		void set_bbox();	// set bounding box
		bool is_in(double xf[2]);
		Ellip();
	private:
	protected:
};
//------------------ Pixel Class ------------------
class Pixel{
	public:
		double Xa[2], Xb[2];
		double Wd[2],Xc[2];
		double xs[4],ys[4];
		double Le[4];
		double rmax,rmin;
		void print();
		void draw(char fn[128],char mode[3]);
		void draw();
		void set_Xa(double xa,double ya);
		void set_Xb(double xb,double yb);
		double dist2pt(double px, double py);
		bool ready;
		bool is_in;
		void setup();
		Pixel();
	private:
	protected:
};
//---------------------------------------------------
int is_cross(Pixel px, Circ cr);
void translate_crs(int icrs);
//---------------------------------------------------
class Solid{
	public:
		int nelp;
		Ellip *els;
		Solid();
		Solid(int n);
		void draw(char fn[128],int ndat);
	private:
	protected:
};
class QPatch{
	public:
		Pixel px;	// pixel data
		QPatch *par;	// parent
		QPatch *chld[4]; // children
		QPatch();	// default constructor
		void print();	// 
		void draw(char fn[128],char mode[3]);
		void draw();
		int lev;
		int icrs;	//  0:exterior/ 1:interior/ 2:inclusive/ 3:boundary pixel
		bool bndr;
		bool intr;
		bool extr;
		void set_lim(double xa[2], double xb[2]);
	private:
	protected:
};
QPatch *new_QPatch(double Xa[2], double Wd[2]);
void gather_leaves(QPatch *qp_par, int *count, QPatch *qp_leaves);
int Qtree(QPatch *qp, Circ cr,int *count);
int Qtree(QPatch *qp, Solid sld,int *count);

class Poly{
	public:
		int np;	
		double *xs,*ys;
		double xg[2];	// centroide
		void draw();
		void draw(char fn[128],char mode[3]);
		void mem_alloc();
		Poly();
		Poly(int n);
		bool is_in(double xf[2]);
		void slide(double ux, double uy);
		void rotate(double xc[2], double th);
		void rotate(int node_num, double th);
		Bbox bbox;	// set bounding box
		void set_bbox();	// set bounding box
		void set_center();
	private:
	protected:
};
int poly_cross(Poly A, Poly B);
int poly_cross(Poly A, Circ B);
int poly_cross(Poly A, Ellip B);
