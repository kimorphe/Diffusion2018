#include<stdio.h>
#include<stdlib.h>
//#include<random>
#include"set2d.h"
//#include"vec2.h"
#ifndef __TCNTRL__
	#define __TCNTRL__
	#include "tcntrl.h"
#endif
#ifndef __GRID__
	#define __GRID__
	#include "grid.h"
#endif

#include "pore.h"
using namespace std;
int main(){
	char fn[128]="pore.dat";
	PoreCells Pcll;
	Pcll.load_cell_data(fn);

	int ng=Pcll.count_grids();
	printf("number of grids=%d\n",ng);
	Grid gd(ng);
	gd.set_grid_params(Pcll.Xa,Pcll.Xb,Pcll.Nx,Pcll.Ny);

	int i,j,iad=0,ID=0;
	for(i=0; i<Pcll.Nx; i++){
	for(j=0; j<Pcll.Ny; j++){
		if(Pcll.grid_type(i,j)==1 || Pcll.grid_loc(i,j)==1){
//		if(Pcll.grid_type(i,j)==1){
		       gd.NDs[iad].id=ID;
		       gd.NDs[iad].iad=iad;
		       iad++;
		};
		ID++;
	}
	}
	gd.connect();

	double xcod,ycod;

	std::mt19937_64 engine(-5);
	std::uniform_real_distribution<double>MT01(0.0,1.0);
	int next,now;
	int nwk=1,Nt=400,inc=1; 
	
	Node **nd0=(Node **)malloc(sizeof(Node*)*nwk);
	FILE *fp=fopen("rw.out","w");
	for(i=0;i<nwk;i++){
		now=int(MT01(engine)*gd.Ng);
		nd0[i]=&gd.NDs[now];
	};

	fprintf(fp,"%d,%d,%d\n",nwk,Nt,inc);
	fprintf(fp,"%le,%le\n",gd.dx[0],gd.dx[1]);
	double x0,y0;
	double tolx=gd.dx[0]*1.001;
	double toly=gd.dx[1]*1.001;
	double *ofx=(double *)malloc(sizeof(double)*nwk);
	double *ofy=(double *)malloc(sizeof(double)*nwk);
	for(j=0;j<Nt;j++){
		printf("step=%d\n",j);
		for(i=0;i<nwk;i++){
			gd.grid_cod(nd0[i]->iad,&x0,&y0);
			next=int(MT01(engine)*4)%4;
			if(nd0[i]->cnct[next]!=-1){
				nd0[i]=nd0[i]->cnds[next];
				gd.grid_cod(nd0[i]->iad,&xcod,&ycod);

				if(xcod-x0>tolx) ofx[i]-=gd.Wd[0];
				if(x0-xcod>tolx) ofx[i]+=gd.Wd[0];
				if(ycod-y0>toly) ofy[i]-=gd.Wd[1];
				if(y0-ycod>toly) ofy[i]+=gd.Wd[1];
				if(j%inc==0) fprintf(fp,"%lf,%lf\n",xcod+ofx[i],ycod+ofy[i]);
			}
		}
	}
	return(0);
};
