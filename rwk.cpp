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

int main(int argc, char *argv[]){
	char fdat[128];	// input  (pore cell data)
	char fout[128]; // output (walker position)
	char fu2b[128]; // output (square mean displacement)
	int nwk;	// number of walkers
	int Nt;	// time steps
	int inc; 	// output time step increment
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
	puts(fdat);
	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fout);
	puts(fout);
	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fu2b);
	puts(fu2b);

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nwk);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&Nt);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&inc);
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
//		printf("(i0,j0)=(%d, %d)\n",i0,j0);
		for(k=0;k<4;k++){
//			printf("k=%d,cnct[k]=%d\n",k,cnct[k]);
			gd.NDs[l].cnct[k]=cnct[k];
			if(cnct[k]==-1) continue;
			i=i0+iofs[k];
			j=j0+jofs[k];
			if(i<0) i+=gd.Nx;
			if(j<0) j+=gd.Ny;
			if(i>=gd.Nx) i-=gd.Nx;
			if(j>=gd.Ny) j-=gd.Ny;
			iad=gd.find(i*gd.Ny+j);
//			printf("ID%d, (i,j)=(%d,%d)\n",i*gd.Ny+j,i,j);
//			printf("ityp=%d\n",Pcll.grid_type_verb(i,j));
			if(iad==-1) puts("ERROR!!");
			gd.NDs[l].cnds[k]=gd.NDs+iad;
			nc++;
		}
	};
	printf("nc(mean)=%lf\n",nc/(float)ng);

//	gd.connect();

	gd.setup_walkers(nwk,-5);	// setup random walkers 
	gd.init_rand(-2);

	FILE *fu=fopen(fu2b,"w");
	for(j=0;j<Nt;j++){
		gd.rwk();
	//	if(j%inc==0) gd.write_wks(fout);
		fprintf(fu,"%lf\n", gd.mean_u2());
	};
	gd.write_wks(fout);
	fclose(fu);

/*
	std::mt19937_64 engine(-5);
	std::uniform_real_distribution<double>MT01(0.0,1.0);
	double xcod,ycod;
	int next,now;
	
	Node **nd0=(Node **)malloc(sizeof(Node*)*nwk);
	fp=fopen(fout,"w");
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
	double xb,yb;
	for(j=0;j<Nt;j++){
		if(j%10==0) printf("step=%d/%d\n",j,Nt);
		xb=0.0; yb=0.0;
		for(i=0;i<nwk;i++){
			gd.grid_cod(nd0[i]->iad,&x0,&y0);
			xcod=x0; ycod=y0;
			next=int(MT01(engine)*4)%4;
			if(nd0[i]->cnct[next]!=-1){
				nd0[i]=nd0[i]->cnds[next];
				gd.grid_cod(nd0[i]->iad,&xcod,&ycod);

				if(xcod-x0>tolx) ofx[i]-=gd.Wd[0];
				if(x0-xcod>tolx) ofx[i]+=gd.Wd[0];
				if(ycod-y0>toly) ofy[i]-=gd.Wd[1];
				if(y0-ycod>toly) ofy[i]+=gd.Wd[1];
				if(j%inc==0) fprintf(fp,"%lf,%lf\n",xcod+ofx[i],ycod+ofy[i]);
				xcod+=ofx[i];
				ycod+=ofy[i];
			}
			xb+=xcod; 
			yb+=ycod;
		}
		xb/=nwk;
		yb/=nwk;
	}
*/
	return(0);
};
