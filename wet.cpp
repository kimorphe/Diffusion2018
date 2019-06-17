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
/*
 	int Lev=9;	// Quad-tree height
	double Wd[2]={1.0,1.0}; // Unit Cell Size
	double Xa[2]={0.0,0.0}; // Unit Cell position (lowerleft vertex)
	char fsld[128]="solid.dat";	// solid phase data file (input)
	char fout[128]="pore.dat";	// pore coverning cell data (output)
	double thE=30.0; // contact angle
	double Sr=0.5;	// degree of saturation
	double T1=1.e0,T2=1.e-06;
	int nstep=200;
*/

	int Lev;	// Quad-tree height
	double Wd[2]; // Unit Cell Size
	double Xa[2]; // Unit Cell position (lowerleft vertex)
	char fsld[128];	// solid phase data file (input)
	char fout[128];	// pore coverning cell data (output)
	double thE; // contact angle
	double Sr;	// degree of saturation
	double T1,T2;
	int nstep;

//	---------------------------------------
	FILE *fp=fopen("wet.inp","r");
	FILE *fo=fopen("wet.erg","w");
	char cbff[128];

	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fsld);
	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fout);

	puts(fsld);
	puts(fout);

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&Lev);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf,%lf\n",Xa,Xa+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf,%lf\n",Wd,Wd+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&thE);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&Sr);
	fgets(cbff,128,fp);
	fscanf(fp,"%le,%le\n",&T1,&T2);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nstep);
	fclose(fp);
//	---------------------------------------

	sld.load(fsld);	// import solid phase data
	sld.bbox.setup(Xa,Wd); // set bounding box


	//Tree4 tr4;
	//tr4.setup(sld.els,sld.nelp,false,Lev,sld.bbox);
	//tr4.draw();

	PoreCells Pcll;
	Pcll.load_gmm(thE);
	Pcll.qp0.refine[0]=true;	// set parameter to refine pore space plus boundary
	Pcll.setup(sld.els,sld.nelp,false,Lev,sld.bbox); // setup pore coverning regular cells 
	Pcll.connect(); // establish connection among pore coverning cells
	Pcll.init(Sr);	// initialize phase distribution
	printf("total energy=%lf\n",Pcll.total_energy());

	Temp_Hist TH(T1,T2,nstep);

	int i=0;
	double dE;
	while(TH.cont_iteration){
		TH.inc_Temp_exp();
		dE=Pcll.MC_stepping(TH);
		if(i%10==0) printf("%d/%d %10.05e %10.05e\n",i,nstep,dE,Pcll.Etot);
		fprintf(fo,"%12.06e %12.06e %12.06e\n",TH.Temp,dE,Pcll.Etot);
		i++;
	};
	Pcll.write_phs();
	Pcll.fwrite_cells(fout);
	Pcll.write_leaves();

	return(0);
};
