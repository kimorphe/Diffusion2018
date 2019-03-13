#define DB 0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>
#include "set2d.h"

using namespace std;

#if DB==0
int main(){
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> RndR(0.0,1.0);

	char fn[128]="log.txt";
	char md[3]="w",mode[3];

	double PI=4.0*atan(1.0);
	int i,np=100;
	double x,y;
	double ra,rb,aspect=0.6;
	double phi;
	double Wd[2];

	Solid sld(np);

	Wd[0]=1.0;
	Wd[1]=1.0;


	for(i=0;i<np;i++){
		x=RndR(mt)*Wd[0];	
		y=RndR(mt)*Wd[1];	
		ra=0.1*Wd[0];
		rb=ra*aspect;
		phi=RndR(mt)*PI;
		//printf("x,y=%lf %lf, phi=%lf[deg]\n",x,y,phi/PI*180.0);
		sld.els[i].set_xc(x,y);
		sld.els[i].phi=phi;
		sld.els[i].set_radi(ra,rb);
		sld.isect[i]=false;
		sld.els[i].set_bbox();
	};

	QPatch qp0;
	double Xa[2];
	Xa[0]=0.0;
	Xa[1]=0.0;
	qp0.set_lim(Xa,Wd);
	int count=0;

	int lev_max=8;
	Qtree(&qp0,sld,&count,lev_max);
	printf("number of leaves =%d\n",count);

	QPatch *qp_leaves=(QPatch *)malloc(sizeof(QPatch)*count);
	count=0;
	gather_leaves(&qp0,&count,qp_leaves);

	char fname[128];
	sprintf(fname,"qtree_in.out");
	sld.draw(fname,100);
	sprintf(mode,"a");
	for(int i=0;i<count;i++){
		if(qp_leaves[i].intr) qp_leaves[i].draw(fname,mode);
	}

	sprintf(fname,"qtree_bnd.out");
	sld.draw(fname,100);
	for(int i=0;i<count;i++){
		if(qp_leaves[i].bndr) qp_leaves[i].draw(fname,mode);
	}

	sprintf(fname,"qtree.out");
	sld.draw(fname,100);
	sprintf(mode,"a");
	for(int i=0;i<count;i++){
		qp_leaves[i].draw(fname,mode);
	}


	return(0);
};
#endif
#if DB==2	// testing Solid class
int main(){

	return(0);	
};
#endif
