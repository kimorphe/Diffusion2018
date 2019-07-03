#include<stdio.h>
#include<stdlib.h>
//#include<random>
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
	char fdat[128];	// input  (pore cell data)
	char fout[128]; // output (walker position)
	char fu2b[128]; // output (square mean displacement)
	int nwk;	// number of walkers
	int Nt;	// time steps
	int nout,inc; 	// output times and step increment
	char cbff[128];
	FILE *fp;

//	----------------------------------
	if(argc==1){
		fp=fopen("rwk.inp","r");
	}else{
		fp=fopen(argv[1],"r");
	};
	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fdat);
	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fout);
	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fu2b);

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nwk);
	fgets(cbff,128,fp);

	fscanf(fp,"%d\n",&Nt);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nout);
	inc=Nt/nout;
	if(inc==0) inc=1;
	printf("nout=%d\n",nout);
	printf("inc=%d\n",inc);
	fclose(fp);
//	----------------------------------


	PoreCells Pcll;
	Pcll.load_cell_data(fdat);
	int ng=Pcll.count_grids();
	printf("number of grids=%d\n",ng);
	Grid gd(ng);
	gd.set_grid_params(Pcll.Xa,Pcll.Xb,Pcll.Nx,Pcll.Ny);

	int i,j,k,iad=0,ID=0;
	for(i=0; i<Pcll.Nx; i++){
	for(j=0; j<Pcll.Ny; j++){
//		if(Pcll.grid_type(i,j)==1 || Pcll.grid_loc(i,j)==1){
		if(Pcll.grid_type(i,j)==1){
		       gd.NDs[iad].id=ID;
		       gd.NDs[iad].iad=iad;
		       iad++;
		};
		ID++;
	}
	}
	printf("ng=%d,iad_final=%d\n",ng,iad);


	int l,i0,j0;
	int iofs[4]={-1, 0, 1, 0};
	int jofs[4]={ 0,-1, 0, 1};
	int cnct[4];
	int nc=0;
	for(l=0;l<ng;l++){
		ID=gd.NDs[l].id;
		i0=ID/gd.Ny;
		j0=ID%gd.Ny;
		Pcll.grid_connect(i0,j0,cnct);
		for(k=0;k<4;k++){
			gd.NDs[l].cnct[k]=cnct[k];
			if(cnct[k]==-1) continue;
			i=i0+iofs[k];
			j=j0+jofs[k];
			if(i<0) i+=gd.Nx;
			if(j<0) j+=gd.Ny;
			if(i>=gd.Nx) i-=gd.Nx;
			if(j>=gd.Ny) j-=gd.Ny;
			iad=gd.find(i*gd.Ny+j);
			if(iad==-1) puts("ERROR!!");
			gd.NDs[l].cnds[k]=gd.NDs+iad;
			nc++;
		}
	};
	printf("nc(mean)=%lf\n",nc/(float)ng);
	//gd.dbg_connectivity();

	gd.setup_walkers(nwk,-5);	// setup random walkers 
	gd.init_rand(-2);

	FILE *fu=fopen(fu2b,"w");
	int iout=0;
	char ftmp[128];
	double uxb,uyb,u2b;
	for(j=0;j<Nt;j++){
		if(j%inc==0){
			sprintf(ftmp,"rwk%d.out",iout);
			gd.write_wks(ftmp,j);
			iout++;
		}
		gd.rwk();
		u2b=gd.mean_u2();
		gd.mean_u(&uxb, &uyb);
		fprintf(fu,"%lf %lf %lf\n", u2b,uxb,uyb);
	};
	gd.write_wks(fout,Nt-1);

	fclose(fu);

	return(0);
};
