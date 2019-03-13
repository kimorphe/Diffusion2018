#define DB 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>
#include "set2d.h"

using namespace std;


#if DB==1
int main(){
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> RndR(0.0,1.0);
	double PI=4.0*atan(1.0);
	int i;
	double x,y;
	double ra,rb,aspect=0.6;
	double phi;
	double Wd[2],Xa[2],Xb[2];

	Xa[0]=-1.0; Xa[1]=-1.0;
	Xb[0]=1.0; Xb[1]=1.0;
	Wd[0]=1.0; Wd[1]=1.0;

	Ellip el1,el2;
	el1.set_xc(0.0,-0.5);
	el1.set_radi(1.0,0.5);
	el1.set_phi(0.0);
	el1.set_bbox();
	printf("A1=%lf\n",el1.area());

	el2.set_xc(0.0,0.5);
	el2.set_radi(1.0,0.5);
	el2.set_phi(0.0);
	el2.set_bbox();
	printf("A2=%lf\n",el2.area());

	bool isect=false;
	int lev_max=6;
	double S=area(el1,el2,lev_max,isect);
	printf("A=%lf, err=%lf\n",S,S-el1.area()-el2.area());

	Bbox bx;	
	if(isect){
		bx=bbox_cross(el1.bbox,el2.bbox);
	}else{
		bx=bbox_union(el1.bbox,el2.bbox);
	}

	QPatch qp0;
	qp0.set_lim(bx.Xa,bx.Xb);
	int count=0;
	Qtree(&qp0,el1,el2,isect,&count,lev_max);
	printf("count=%d\n",count);
	QPatch *qp_leaves=(QPatch *)malloc(sizeof(QPatch)*count);
	count=0;
	gather_leaves(&qp0,&count,qp_leaves);

	char mode[3]="w",fname[128];
	sprintf(fname,"qtree_in.out");
	el1.draw(fname,100,mode);
	sprintf(mode,"a");
	el2.draw(fname,100,mode);
	for(int i=0;i<count;i++){
		if(qp_leaves[i].intr) qp_leaves[i].draw(fname,mode);
	}

	sprintf(fname,"qtree_bnd.out");
	sprintf(mode,"w");
	el1.draw(fname,100,mode);
	sprintf(mode,"a");
	el2.draw(fname,100,mode);
	for(int i=0;i<count;i++){
		if(qp_leaves[i].bndr) qp_leaves[i].draw(fname,mode);
	}

	clear_Qtree(&qp0);

	return(0);
};
#endif
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

	int q,p=int(RndR(mt)*np);
	bool pq;
	Ellip ep=sld.els[p];
	ep.draw(50);
	ep.bbox.draw();
	for(q=0;q<np;q++){
		pq=bbox_cross(ep,sld.els[q]);
		if(pq){
		       	//sld.els[q].draw(50);
		       	sld.els[q].bbox.draw();
		}
	}	


	double Xa[2];
	Xa[0]=0.0;
	Xa[1]=0.0;
	QPatch qp0;
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
