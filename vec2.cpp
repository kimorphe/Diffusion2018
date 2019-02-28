#include <stdio.h>
#include <math.h>

//double dist(double *x, double *y);
//void show_msg(char *fname);

class Vec2{
	public:
	double x[2];
	void set(double,double);	
	void set(double*);	
	Vec2 times(double);	
	Vec2 div(double);	
	double len();
	void print();
	private:
};
double iprod(Vec2 a,Vec2 b);
double iprod(double a[2], double b[2]);
Vec2 vsum(Vec2 a, Vec2 b);
Vec2 vdiff(Vec2 a, Vec2 b);
//--------- Vec2 Class Methods ----------
void Vec2::set(double X, double Y){
	x[0]=X; x[1]=Y;
};
void Vec2::set(double X[2]){
	x[0]=X[0]; x[1]=X[1];
};
Vec2 Vec2::times(double s){
	Vec2 v;
	v.set(x[0]*s, x[1]*s);
	return(v);
};
Vec2 Vec2::div(double s){
	Vec2 v;
	v.set(x[0]/s, x[1]/s);
	return(v);
};
double Vec2::len(){
	return(sqrt(x[0]*x[0]+x[1]*x[1]));
};
void Vec2::print(){
	printf("(%lf, %lf)\n",x[0],x[1]);
};

//---- Functions invovling Vec2 Class ------- 
double iprod(Vec2 a,Vec2 b){
	return(a.x[0]*b.x[0]+a.x[1]*b.x[1]);
};
double iprod(double a[2], double b[2]){
	return(a[0]*b[0]+a[1]*b[1]);
};
Vec2 vsum(Vec2 a, Vec2 b){
	Vec2 v;
	v.set(a.x[0]+b.x[0], a.x[1]+b.x[1]);
	return(v);
};
Vec2 vdiff(Vec2 a, Vec2 b){
	Vec2 v;
	v.set(a.x[0]-b.x[0], a.x[1]-b.x[1]);
	return(v);
};
