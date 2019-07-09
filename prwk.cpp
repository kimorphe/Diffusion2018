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

void show_msg(char fnm[128]){
	printf("Error !! Can't open file: %s\n",fnm);
	exit(-1);
};

int main(int argc, char *argv[]){
	char fdat[128];	// input  (pore cell data)
	char fsld[128]; // solid phase data
	char fout[128]; // output (walker position)
	char fu2b[128]; // output (square mean displacement)
	int nwk;	// number of walkers
	int Nt;	// time steps
	int nout,inc; 	// output times and step increment
	char cbff[128];
	FILE *fp;
	int ngap;	// 0:close ,1:open thorat model

	Solid sld;

//	----------------------------------
	char fnm[128];
	if(argc==1){
		sprintf(fnm,"%s","rwk.inp");
		fp=fopen(fnm,"r");
	}else{
		sprintf(fnm,"%s",argv[1]);
		fp=fopen(argv[1],"r");
	};
	if(fp==NULL) show_msg(fnm);
	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fdat);
	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fsld);
	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fout);
	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fu2b);

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&ngap);
	printf("ngap=%d\n",ngap);
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
	sld.load(fsld);	// import solid phase data


	PoreCells Pcll;
	Pcll.load_cell_data(fdat);

	sld.bbox.setup(Pcll.Xa,Pcll.Wd); // set bounding box
	Pcll.connect4();


	Pcll.setup_walkers(nwk,-5);
	//gd.init_rand(-2);
	FILE *fu=fopen(fu2b,"w");
	int j,iout=0;
	char ftmp[128];
	double uxb,uyb,u2b;
	for(j=0;j<Nt;j++){
		if(j%inc==0){
			sprintf(ftmp,"rwk%d.out",iout);
			Pcll.write_wks(ftmp,j);
			iout++;
		}
		Pcll.rwk(-1);
		u2b=Pcll.mean_u2();
		Pcll.mean_u(&uxb, &uyb);
		fprintf(fu,"%lf %lf %lf\n", u2b,uxb,uyb);
	};
	Pcll.write_wks(fout,Nt-1);

	fclose(fu);

	return(0);
};

