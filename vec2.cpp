#define DEBUG 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vec2.h"

//--------- Vec2 Constructor  ----------
Vec2::Vec2(){};
Vec2::Vec2(double xf[2]){
	x[0]=xf[0];
	x[1]=xf[1];
};
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

void Vec2::rotate(double th){
	double xx,yy;
	double cost=cos(th),sint=sin(th);
	xx= x[0]*cost-x[1]*sint;
	yy= x[0]*sint+x[1]*cost;

	x[0]=xx;
	x[1]=yy;
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
double vdist(Vec2 a, Vec2 b){
	Vec2 c;
	c=vdiff(b,a);
	return(c.len());
};
double distP2L(
	double xpt[2],
	double x1[2],
	double x2[2],
	double *rmin,
	double *rmax
){
	Vec2 y1(x1);
	Vec2 y2(x2);
	Vec2 xp(xpt);
	Vec2 tx,nx;

	*rmax=0.0;

	double r1,r2,rp,r12;
	
	r1=vdist(y1,xp);
	r2=vdist(y2,xp);

	tx=vdiff(y2,y1);
	r12=tx.len();
	tx.div(r12);

	nx.set(-tx.x[1],tx.x[0]);

	*rmin=r1;
	*rmax=r1;

	if(*rmin > r2) *rmin=r2;
	if(*rmax < r2) *rmax=r2;

	double Xt=iprod(tx,vdiff(xp,x1));

	if(Xt < 0.0) return(*rmin);
	if(Xt > r12) return(*rmin);

	double Xn=iprod(nx,vdiff(xp,x1));
	
	*rmin=fabs(Xn);
	return(*rmin);
};
#if DEBUG == 1
int main(){
	Vec2 x1,x2,xp;
	x1.set(1.0,0.0);
	x2.set(2.0,0.0);

	double xx,yy;
	double rmin,rmax;
	for(int i=0;i<=30;i++){
		xx=i*0.1;
		yy=i*0.1-1.0;;
		xp.set(0.0,yy);
		distP2L(xp.x, x1.x, x2.x,&rmin,&rmax);
		printf("%lf %lf %lf %lf\n",xx,yy,rmin,rmax);
	};

	return(0);
};
#endif
