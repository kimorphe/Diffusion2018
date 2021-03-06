#include<stdio.h>
#include<stdlib.h>
#include<time.h>
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
int main(int argc, char *argv[]){
	Solid sld;

	int Lev;	// Quad-tree height
	double Wd[2]; // Unit Cell Size
	double Xa[2]; // Unit Cell position (lowerleft vertex)
	char fsld[128];	// solid phase data file (input)
	char fout[128];	// pore coverning cell data (output)
	char foutc[128];	// pore coverning cell data (output)
	double thE; // contact angle
	double Sr;	// degree of saturation
	double T1,T2;
	int nstep;
	int seed;
	int rstat;	// restart (0:No, 1: yes)
	char f_rstat[128]; // file name (import data for restar MC)

//	---------------------------------------
	FILE *fp;
	char fnm[128];
	if(argc >1){
		sprintf(fnm,"%s",argv[1]);
	}else{
		sprintf(fnm,"%s","wet.inp");
	}
	fp=fopen(fnm,"r");
	if(fp==NULL){
		printf("Error !! Can't open file:%s\n",fnm);
		exit(-1);
	};
	FILE *fo=fopen("wet.erg","w");
	char cbff[128];
	double tmp;

	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fsld);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&rstat);
	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",f_rstat);
	printf("%s\n",f_rstat);

	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fout);
	fscanf(fp,"%s\n",foutc);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&seed);
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
	double rstep;
	fscanf(fp,"%lf\n",&rstep);
	fclose(fp);
//	---------------------------------------

	sld.load(fsld);	// import solid phase data
	sld.bbox.setup(Xa,Wd); // set bounding box

	clock_t tstart,tend; 

	tstart=clock();

	char ftmp[128]="pore_initial.dat";
	double Etot0,Etot;
	PoreCells Pcll,Pcllc;
	Pcll.load_gmm(thE);	// set contact angle thE
	Pcll.ngap=2;		// set inter-particle gap (1:close,2:open) 
	Pcll.qp0.refine[0]=true;// set parameter to refine pore space plus boundary

	int Lev_Exact=10;
	Pcll.isetup(sld.els,sld.nelp,false,Lev,Lev_Exact,sld.bbox); // setup pore coverning regular cells 
	//Pcll.fwrite_cells("log_cell.dat");
	//Pcll.setup(sld.els,sld.nelp,false,Lev,sld.bbox); // setup pore coverning regular cells 

	Pcllc.load_gmm(thE);
	Pcllc.ngap=1;	// close thorat
	Pcllc.qp0.refine[0]=true;
	Pcllc.Sr=Sr;
	//Pcllc.setup(sld.els,sld.nelp,false,Lev,sld.bbox); // setup pore coverning regular cells 
	Pcllc.isetup(sld.els,sld.nelp,false,Lev,Lev_Exact,sld.bbox); // setup pore coverning regular cells 
	Pcllc.connect(); // establish connection among pore coverning cells


	Pcll.connect(); // establish connection among pore coverning cells
	if(rstat==0){ //fresh start
		Pcll.init(Sr);	// initialize phase distribution
	}else{	// restart
		//Pcll.fwrite_cells("error0.dat");
		PoreCells pc0;
		pc0.load_cell_data(f_rstat);
		Pcll.Sr=pc0.Sr;
		puts("Copying data...");
		refine_PoreCell_data(&pc0,&Pcll);
		puts("Done");
		Pcll.setup_phs_list();
	};
	Pcll.fwrite_cells(ftmp);

	nstep=Pcll.n_water*rstep;
	printf("Total MC steps =%d\n",nstep);
	Temp_Hist TH(T1,T2,nstep);

	int i,j,jmax=100;
	int nswap,nswap_sum;
	double dE,Evar;
	double alph=0.05;


	Etot0=Pcll.total_energy();
	printf("Total energy=%lf\n",Etot0);
	for(j=0;j<jmax;j++){
		i=0;
		Evar=0.0;
		nswap_sum=0;
		//printf("j=%d\n",j);
		while(TH.cont_iteration){
			TH.inc_Temp_exp();
			dE=Pcll.MC_stepping(TH,&nswap,seed);
			Evar+=dE*dE;
			fprintf(fo,"%12.06e %12.06e %12.06e\n",TH.Temp,dE,Pcll.Etot);
			nswap_sum+=nswap;
			i++;
		};
		Evar=sqrt(Evar)/nstep;
		printf("%d: Estd=%le, nswap=%d \n",j,Evar,nswap_sum);
		if(Evar<1.e-05) break;
		TH.renew(alph);
	};


	Etot=Pcll.total_energy();
	printf("Total energy=%lf\n",Etot);
	Pcll.write_phs();
	Pcll.fwrite_cells(fout);
	Pcll.write_leaves();

	copy_PoreCell_data(&Pcll,&Pcllc);
	//sprintf(fout,"pore_c.dat");
	Pcllc.fwrite_cells(foutc);
	tend=clock();

	printf("Elapsed time= %lf[sec]\n ",double(tend-tstart)/CLOCKS_PER_SEC);

	return(0);
};
