#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>
#ifndef __TCNTRL__
	#define __TCNTRL__
	#include "tcntrl.h"
#endif
using namespace std;

//--------------Bounding Box  Class--------------------
class Bbox{
	public:
	double Xa[2];
	double Xb[2];
	double Wd[2];
	void set_Xa(double x, double y);
	void set_Xb(double x, double y);
	void set_Wd();
	void setup(double xa[2], double xb[2]);
	void draw();
	void draw(char fn[128],char mode[3]);
	void slide(double ux, double uy);
	double area();
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
		void slide(double ux,double uy, Bbox unit_cell );
		void slide(double ux,double uy);
		Bbox bbox;	// bounding box
		void set_bbox();	// set bounding box
		bool is_in(double xf[2]);
		Ellip();
		double area();
	private:
	protected:
};
bool bbox_cross(Ellip el1, Ellip el2);
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
		void draw(FILE *fp);
		void write(FILE *fp);
		void draw();
		void set_Xa(double xa,double ya);
		void set_Xb(double xb,double yb);
		double dist2pt(double px, double py);
		bool ready;
		bool is_in;
		void setup();
		double area();
		Pixel();
	private:
	protected:
};
//---------------------------------------------------
bool bbox_cross(Bbox b1, Pixel px);
int is_cross(Pixel px, Circ cr);
void translate_crs(int icrs);
//---------------------------------------------------
class Solid{
	public:
		int nelp;
		bool *isect;	// true for set intersection, faulse for union 
		Ellip *els;
		Solid();
		Solid(int n);
		Solid(int n, double Wd[2]);
		double A0;
		void draw(char fn[128],int ndat);
		Bbox bbox;
		double perturb(int p, double ux, double uy, double dphi);
		double perturb_periodic(int p, double ux, double uy, double dphi);
		double MC(Temp_Hist TH);
		void init_rand(int seed);
		double area(int Lev_Max);
		void load(char fn[128]);
		void write(char fn[128]);
		std::mt19937_64 mt;	// random number generator
		std::uniform_real_distribution<> Urnd;
		std::normal_distribution<> Grnd;
		void set_domain(double Xa[2], double Wd[2]);
		double s_over;	// overlapped area ratio
		double poro;
	private:
	protected:
};
class QPatch{
	public:
		Pixel px;	// pixel data
		QPatch *par;	// parent
		QPatch *chld[4]; // children
		QPatch();	// default constructor
		~QPatch();	// destructor
		void print();	// 
		void draw(char fn[128],char mode[3]);
		void draw(FILE *fp);
		void write(FILE *fp);
		void draw();
		int lev;
		int icrs;	//  0:exterior/ 1:interior/ 2:inclusive/ 3:boundary pixel
		bool refine[4]; // Patch is/isn't to be refined if refine[icrs] = true/false 
		bool bndr;
		bool intr;
		bool extr;
		void set_lim(double xa[2], double xb[2]);
		int isin();	// 0:interior, 1:boundary, 2:exterior, -1:error
		void show_prms();
	private:
	protected:
};
QPatch *new_QPatch(double Xa[2], double Wd[2]);
void new_QPatch(double Xa[2], double Wd[2], QPatch **qp_new);
void gather_leaves(QPatch *qp_par, int *count, QPatch *qp_leaves);
int Qtree(QPatch *qp, Circ cr,int *count);
int Qtree(QPatch *qp, Solid sld,int *count, int lev_max);
int Qtree(QPatch *qp, Ellip el1, Ellip el2, bool isect, int *count, int lev_max);
int Qtree(QPatch *qp, Ellip *els, int nelp, bool isect, int *count, int lev_max);
int Qtree(QPatch *qp, Ellip *els, int nelp, bool isect, int *count, int lev_max, Bbox unit_cell);
int QtreeFind(QPatch *qp, double xf[2]);
double area(Ellip el1, Ellip el2, int lev_max, bool isect);
void clear_Qtree(QPatch *qp);
void clear_Qtree2(QPatch *qp);
void indx4PrdBC(		// retrun indices to apply periodic B.C.
	int *i1, int *i2, 	// 1st index (start,end)
	int *j1, int *j2, 	// 2nd index (start,end)
	Bbox bx, 		// bounding box of geoemetric object
	Bbox B0			// Unit cell
);
class Tree4{
	public:
		Tree4();
		QPatch qp0;	// Root node
		bool isect;	// true:intersection, false:union
		bool ready;	// Has tree been build already?
		int lev_max;	// deepest level
		int n_leaves;	// number of leaves
		void setup(Ellip elp1, Ellip elp2, bool set_opr, int Lev_Max);
		void setup(Ellip *els, int nelp,bool set_opr, int Lev_Max);
		void setup(Solid sld, int Lev_Max);
		void setup(Ellip *els, int nelp,bool set_opr, int Lev_Max,Bbox bx);
		QPatch *leaves;
		void draw();
		void write();
		void clean();
		double area();
		int nint,next,nbnd;
		void count();
		void set_grid_params();
		void write_grid_params(FILE *fp);
		int Nx,Ny;
		double dx[2],Xa[2],Xb[2],Wd[2];
		int grid_type(int id, int jd, int cnct[4]);
	private:
};
class Poly{
	public:
		int np;	
		double *xs,*ys;
		double xg[2];	// centroide
		void draw();
		void draw(char fn[128],char mode[3]);
		void mem_alloc();
		void mem_free();
		bool alocd;
		Poly();
		Poly(int n);
		~Poly();
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
