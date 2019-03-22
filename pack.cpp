#define DB 0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>
#include "set2d.h"
#ifndef __TCNTRL__
	#define _TCNTRL
	#include "tcntrl.h"
#endif

using namespace std;

#if DB==3	// Debugging memory leak 
// -->  memory blocks for Poly::xs,ys are not freed 
int main(){
	int np=200;
	double Wd[2]={1.0,1.0};
	char fn[128]="geom0.dat";

	Solid sld(np,Wd);
	sld.draw(fn,50);

	Tree4 tr;
	for(int i=0;i<1;i++){
		printf("i_try=%d\n",i);
		tr.setup(sld.els,sld.nelp,false,10,sld.bbox);
//		tr.setup(sld,10);
		tr.draw();
		tr.clean();
	};

	int lev_max=9;
	printf("S=%lf\n",sld.area(lev_max));

/*
	Tree4 tr4;
	tr4.setup(sld,9);
	tr4.draw();
	tr4.clean();
*/
};
#endif

#if DB==0
int main(){
	//std::uniform_real_distribution<double> RndR;
	//RndR=uniform_real_distribution<double>(0.0,1.0);
	Temp_Hist TH;
	char fnt[128]="temp_hist.inp";
	TH.load(fnt);

	int np=500;
	double Wd[2]={1.0,1.0};
	char fn[128]="geom0.dat";
	double dE_tot=0.0,dE;
	Solid sld(np,Wd);
	sld.draw(fn,50);

	sld.area(9);
	while(TH.cont_iteration){
		dE=sld.MC(TH);
		dE_tot+=dE;
		printf("dE=%lf\n",dE);
		TH.inc_Temp();
	};
	printf("dE_tot=%lf\n",dE_tot);
	sprintf(fn,"geom1.dat");
	sld.draw(fn,50);
	sld.area(9);

	Tree4 tr4;
	tr4.setup(sld.els,sld.nelp,false,9,sld.bbox);
	tr4.draw();
	tr4.clean();
	return(0);
};
#endif


#if DB==1	
// checking area(el1,el2,levmax,isect)
// Introducing Tree4 class with which quad-tree based quadrature can be done 
int main(){

	Ellip el1,el2;
	el1.set_xc(0.0,-0.0);
	el1.set_radi(1.0,0.5);
	el1.set_phi(0.0);
	el1.set_bbox();
	printf("A1=%lf\n",el1.area());

	el2.set_xc(0.0,0.5);
	el2.set_radi(1.0,0.5);
	el2.set_phi(0.0);
	el2.set_bbox();
	printf("A2=%lf\n",el2.area());

	bool isect=true;
	int lev_max=6;
	double S=area(el1,el2,lev_max,isect);
	printf("A=%lf, err=%lf\n",S,S-el1.area()-el2.area());

	Tree4 tr4;
	tr4.setup(el1,el2,isect,lev_max);
	tr4.draw();
	printf("A=%lf\n",tr4.area());
	tr4.clean();

	char fname[128]="geom.dat",md[3]="w";	
	el1.draw(fname,100,md);
	sprintf(md,"a");	
	el2.draw(fname,100,md);


	return(0);
};
#endif
#if DB==2
int main(){
	std::random_device rd;
	//std::mt19937 mt(rd());
	std::mt19937_64 mt(-1);
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

	sprintf(fn,"geom.dat");
	sld.draw(fn,50);

	int q,p=int(RndR(mt)*np);
	bool pq;
	Ellip elp=sld.els[p],elq;
	//elp.draw(50);
	//ep.bbox.draw();
	double dE=0.0,S=0.0;
	Tree4 tr4;
	for(q=0;q<np;q++){
		if(q==p) continue;
		pq=bbox_cross(elp,sld.els[q]);
		if(pq){
			elq=sld.els[q];
		       	//elq.draw(50);
		       	//sld.els[q].bbox.draw();
			tr4.setup(elp,elq,true,6);
			S=tr4.area();
			tr4.clean();
			printf("S=%lf\n",S);
			dE+=S;
		}
	}	
	printf("dE=%lf\n",dE);

	double Xa[2],Xb[2];
	Xa[0]=0.0; 
	Xa[1]=0.0;
	Xb[0]=Xa[0]+Wd[0];
	Xb[1]=Xa[1]+Wd[1];
	sld.bbox.setup(Xa,Xb);

	Ellip elr=elp;
	printf("xc=%lf %lf -->",elr.xc[0],elr.xc[1]);
	double ux,uy;
	ux=RndR(mt)-0.45;
	uy=RndR(mt)-0.15;
	elr.slide(ux,uy,sld.bbox);
	printf("xc=%lf %lf\n",elr.xc[0],elr.xc[1]);

	dE=0.0; S=0.0;
	for(q=0;q<np;q++){
		if(q==p) continue;
		pq=bbox_cross(elr,sld.els[q]);
		if(pq){
			elq=sld.els[q];
			tr4.setup(elr,elq,true,6);
			S=tr4.area();
			tr4.clean();
			printf("S=%lf\n",S);
			dE+=S;
		}
	}	
	printf("dE=%lf\n",dE);

	sld.perturb(p,ux,uy,0.0);

	return(0);
};
#endif
