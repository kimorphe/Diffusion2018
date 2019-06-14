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
	Solid sld;

	int Lev=9;	// Quad-tree height
	double Wd[2]={1.0,1.0}; // Unit Cell Size
	double Xa[2]={0.0,0.0}; // Unit Cell position (lowerleft vertex)
	char fn[128]="solid.dat";	// solid phase data file
	char fdat[128]="pore.dat";

	sld.load(fn);	// import solid phase data
	sld.bbox.setup(Xa,Wd); // set bounding box

	PoreCells Pcll;
	Pcll.load_gmm(30.0);
	Pcll.qp0.refine[0]=true;	// set parameter to refine pore space plus boundary
	Pcll.setup(sld.els,sld.nelp,false,Lev,sld.bbox); // setup pore coverning regular cells 
	//Pcll.connect0(); // establish connection among pore coverning cells
	Pcll.connect(); // establish connection among pore coverning cells
	double Sr=0.5;	// degree of saturation
	Pcll.init(Sr);	// initialize phase distribution
	printf("total energy=%lf\n",Pcll.total_energy());


	double T1=1.e0,T2=1.e-06;
	int nstep=100;
	Temp_Hist TH(T1,T2,nstep);

	double dE;
	while(TH.cont_iteration){
		dE=Pcll.MC_stepping(TH);
		printf("%lf %lf %lf\n",TH.Temp,dE,Pcll.Etot);
		TH.inc_Temp_exp();
	};
	Pcll.write_phs();
	Pcll.fwrite_cells(fdat);


/*
	int ng=Pcll.count_grids();
	printf("number of grids=%d\n",ng);
	Grid gd(ng);
	gd.set_grid_params(Pcll.Xa,Pcll.Xb,Pcll.Nx,Pcll.Ny);

	int i,j,iad=0,ID=0;
	for(i=0; i<Pcll.Nx; i++){
	for(j=0; j<Pcll.Ny; j++){
		if(Pcll.grid_type(i,j)==1){
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
	printf("tols=%lf %lf\n",tolx,toly);
	for(j=0;j<Nt;j++){
		printf("step=%d\n",j);
		for(i=0;i<nwk;i++){
			gd.grid_cod(nd0[i]->iad,&x0,&y0);

			next=int(MT01(engine)*4)%4;

			printf("nc=%d\n",nd0[i]->nc);
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
	//Pcll.draw();
*/
	return(0);
};
