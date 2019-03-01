#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "set2d.h"

int main(){
	
	double rmin,x,y,dx,dy;
	Pixel px;
	px.set_Xa(0.0,0.0);
	px.set_Xb(2.0,1.0);

	double x1=2.0,x2=2.0;
	double y1= -1.0,y2=2.0;
	int i,ndivx=1,ndivy=100;
	dx=0.0; dy=0.0;
	dx=(x2-x1)/ndivx;
	dy=(y2-y1)/ndivy;
	//for(i=0;i<=ndivx;i++){
	for(i=0;i<=ndivy;i++){
		x=x1+dx*i;
		y=y1+dy*i;
		rmin=px.dist2pt(x,y);
		printf("%lf %lf %lf\n",x,y,rmin);
	};

	Circ cr;
	cr.xc[0]=-1.00;
	cr.xc[1]=-0.5;
	cr.radi=1.0;

	char fnc[128]="circ.dat";
	char fnp[128]="pix.dat";
	char mode[3]="w";

	cr.draw(fnc,100,mode);
	px.draw(fnp,mode);
	int icrs;

	icrs=is_cross(px,cr);
	printf("icrs=%d\n",icrs);


	return(0);
};
