//#define DB 0
//#define DB 4
#define DB -1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>
#include "set2d.h"
#ifndef __TCNTRL__
	#define __TCNTRL__
	#include "tcntrl.h"
#endif

#ifndef __GRIDL__
	#define __GRID__
	#include "grid.h"
#endif

using namespace std;
Node::Node(){
	id=-1;
	for(int i=0;i<4;i++){
	       	cnct[i]=-1;
		cnds[i]=NULL;
	}
};
int Grid::find(int k){
	if(k<NDs[0].id) return(-1);
	if(k>NDs[Ng-1].id) return(-1);
	int i1=0;
	int i2=Ng-1;
	int im;

	if(NDs[i1].id==k) return(i1);
	if(NDs[i2].id==k) return(i2);
	while(i2-i1>1){
		im=(i1+i2)/2;
		if(NDs[im].id==k) return(im);
		if(NDs[im].id>k){
			i2=im;
		}else{
			i1=im;
		}
	}

	if(NDs[i1].id==k) return(i1);
	if(NDs[i2].id==k) return(i2);
	return(-1);
};

void Grid::setup(Tree4 tr4){
	tr4.set_grid_params();

	Nx=tr4.Nx;
	Ny=tr4.Ny;
	for(int i=0;i<2;i++){
		Xa[i]=tr4.Xa[i];
		Xb[i]=tr4.Xb[i];
		Wd[i]=tr4.Wd[i];
		dx[i]=tr4.dx[i];
	}

	int cnct[4];
	int ngb=0;
	int m,k;
	for(k=0; k<tr4.Nx; k++){
	for(m=0; m<tr4.Ny; m++){
		if(tr4.grid_type(k,m,cnct)==1) ngb++;
	}
	}
	printf("ngb=%d\n",ngb);
	Ng=ngb;
	Grid::mem_alloc();

	ngb=0;
	int isum=0,l,nc;
	for(k=0; k<tr4.Nx; k++){
	for(m=0; m<tr4.Ny; m++){
		if(tr4.grid_type(k,m,cnct)==1){
			NDs[ngb].id=isum;
			NDs[ngb].iad=ngb;
			nc=0;
			for(l=0;l<4;l++){
			       NDs[ngb].cnct[l]=cnct[l];
			       nc+=cnct[l];
			};
			NDs[ngb].nc=nc;
			ngb++;
		};
		isum++;
	}
	}

	int iofst[4]={-1, 0, 1, 0};
	int jofst[4]={ 0,-1, 0, 1};
	int ix0,jy0,ix,jy,adr,id,ig;
	for(k=0;k<Ng;k++){
		ix0=NDs[k].id/Ny;
		jy0=NDs[k].id%Ny;
		nc=0;
		for(l=0;l<4;l++){
			if(NDs[k].cnct[l]==1){
				ix=ix0+iofst[l];		
				jy=jy0+jofst[l];		
				while(ix<0) ix+=Nx;
				while(jy<0) jy+=Ny;
				while(ix>=Nx) ix-=Nx;
				while(jy>=Ny) jy-=Ny;
				id=ix*Ny+jy;
				adr=Grid::find(id);
				NDs[k].cnds[nc]=(NDs+adr);
				nc++;
			};
		};
	};

};
void Grid::connect(){
	int iad,id,jad;
	int i,j,k,l,ix,iy;
	int ic,ID;
	int iofst[4]={-1, 0, 1, 0};
	int jofst[4]={ 0,-1, 0, 1};
	for(iad=0;iad<Ng;iad++){
		ID=NDs[iad].id;
		i=ID/Ny;
		j=ID%Ny;
		NDs[iad].nc=0;
		ic=0;
		for(ic=0;ic<4;ic++){
			ix=i+iofst[ic];
			iy=j+jofst[ic];
			if(ix<0) ix+=Nx;
			if(ix>=Nx) ix-=Nx;
			if(iy<0) iy+=Ny;
			if(iy>=Ny) iy-=Ny;
			id=ix*Ny+iy;
			jad=Grid::find(id);

			NDs[iad].cnct[ic]=jad;
			if(jad!=-1){
				NDs[iad].cnds[ic]=NDs+jad;
				NDs[iad].nc++;
			};
		}
	};
};
void Grid::l2cod(int l,double *x,double *y){
	int i,j;
	Grid::l2ij(l,&i,&j);
	Grid::indx2cod(i,j,x,y);
};
void Grid::indx2cod(int i,int j,double *x,double *y){
	*x=Xa[0]+dx[0]*i;
	*y=Xa[1]+dx[1]*j;
};
void Grid::l2ij(int l, int *i, int *j){
	*i=l/Ny;
	*j=l%Ny;
};
Grid::Grid(){
	Ng=1;
	ready=false;
};
Grid::Grid(int n){
	Ng=n;
	mem_alloc();
};
void Grid::set_grid_params(double xa[2], double xb[2], int nx, int ny){

	for(int i=0;i<2;i++){
		Xa[i]=xa[i];
		Xb[i]=xb[i];
		Wd[i]=Xb[i]-Xa[i];
	}
	Nx=nx; 
	Ny=ny;
	dx[0]=Wd[0]/Nx;
	dx[1]=Wd[1]/Ny;
	tolx=dx[0]*1.001;
	toly=dx[1]*1.001;
};
void Grid::mem_alloc(){
	NDs=(Node *)malloc(sizeof(Node)*Ng);
	ready=true;
};
void Grid::grid_cod(int inod,double *xcod, double *ycod){
	int id=NDs[inod].id;
	int i,j;
	Grid::l2ij(id,&i,&j);
	Grid::indx2cod(i,j,xcod,ycod);
};
void Grid::setup_walkers(int n, int seed){
	nwk=n;
	wks=(Walker *)malloc(sizeof(Walker)*nwk);

	std::mt19937_64 engine(seed);
	std::uniform_real_distribution<double>MT01(0.0,1.0);
	int i,iad;
	double x0,y0;
	for(i=0;i<nwk;i++){
		iad=int(MT01(engine)*Ng)%Ng;
		wks[i].nd0=&NDs[iad];
		wks[i].ofx=0;
		wks[i].ofy=0;
		grid_cod(wks[i].nd0->iad,&x0,&y0);
		wks[i].x0=x0;
		wks[i].y0=y0;
		wks[i].xn=x0;
		wks[i].yn=y0;
	};
};
void Grid::init_rand(int seed){
	engine=std::mt19937_64(seed);
	irnd3=std::uniform_int_distribution<int>(0,3);
};
void Grid::rwk(){
	Walker wk;
	double xcod,ycod;
	int iwk,next;
	for(iwk=0;iwk<nwk;iwk++){
		wk=wks[iwk];
		next=irnd3(engine);	
		if(wk.nd0->cnct[next]!=-1){
			wk.nd0=wk.nd0->cnds[next];
			grid_cod(wk.nd0->iad,&xcod,&ycod);

			if(xcod-wk.xn > tolx) wk.ofx--;
			if(wk.xn-xcod > tolx) wk.ofx++;
			if(ycod-wk.yn > toly) wk.ofy--;
			if(wk.yn-ycod > toly) wk.ofy++;
			wk.xn=xcod;
			wk.yn=ycod;
			wks[iwk]=wk;
		}
	};
};
void Grid::write_wks(char fname[128]){
	static FILE *fp=fopen(fname,"w");
	Walker wk;
	for(int iwk=0; iwk<nwk; iwk++){
		wk=wks[iwk];
		fprintf(fp,"%lf %lf\n",wk.xn+wk.ofx*Wd[0], wk.yn+wk.ofy*Wd[1]);
	}
};
double Grid::mean_u2(){
	double u2=0.0,dux,duy;
	Walker wk;
	for(int iwk=0; iwk<nwk; iwk++){
		wk=wks[iwk];
		dux=(wk.xn+wk.ofx*Wd[0]-wk.x0);
		duy=(wk.yn+wk.ofy*Wd[1]-wk.y0);
		u2+=(dux*dux+duy*duy);
	};
	return(u2/nwk);
};
//----------------------------------------------------------
