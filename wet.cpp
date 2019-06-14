#include<stdio.h>
#include<stdlib.h>
#include"set2d.h"
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
	Pcll.connect(); // establish connection among pore coverning cells
	double Sr=0.5;	// degree of saturation
	Pcll.init(Sr);	// initialize phase distribution
	printf("total energy=%lf\n",Pcll.total_energy());


	double T1=1.e0,T2=1.e-06;
	int nstep=200;
	Temp_Hist TH(T1,T2,nstep);

	double dE;
	while(TH.cont_iteration){
		dE=Pcll.MC_stepping(TH);
		printf("%lf %lf %lf\n",TH.Temp,dE,Pcll.Etot);
		TH.inc_Temp_exp();
	};
	Pcll.write_phs();
	Pcll.fwrite_cells(fdat);

	return(0);
};
