#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//------------------------Circ Class--------------------
class Circ{
	public:
		double xc[2];	// center 
		double radi;	// radius
		void draw(char fn[128],int npnt, char mode[3]);
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
		void set_lim(double xa[2], double xb[2]);
	private:
	protected:
};
QPatch *new_QPatch(double Xa[2], double Wd[2]);
void gather_leaves(QPatch *qp_par, int *count, QPatch *qp_leaves);
int Qtree(QPatch *qp, Circ cr,int *count);
