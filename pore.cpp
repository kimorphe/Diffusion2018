#include<stdio.h>
#include<stdlib.h>
#include<random>
#include"set2d.h"
#include"vec2.h"
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
//---------------------------------------------------------------------------------
void Material::load(
		double th // contact angle in deg
){
	double PI=4.0*atan(1.0);
	thE=th/180.*PI;	// contact angle in rad

//		0:air, 1:water, 2: solid

	double g10=70.0;	// air-water
	double g20=420.0;	// air-solid
	double g21=g20-g10*cos(thE); // water-solid

	gmm[0][0]=0.0; 
	gmm[1][0]=g10; gmm[1][1]=0.0;
	gmm[2][0]=g20; gmm[2][1]=g21; gmm[2][2]=0.0;

	gmm[0][1]=gmm[1][0];
	gmm[0][2]=gmm[2][0];
	gmm[1][2]=gmm[2][1];

	gmm_max=g10;
	if(gmm_max< g20) gmm_max=g20;

};
void Material::normalize(){
	int i,j;
	for(i=0; i<3; i++){
	for(j=0; j<3; j++){
		gmm[i][j]/=gmm_max;
	}
	}
};
void Material::print_gmm(){
	int i,j;
	printf("-------- Interfacial Energy gmm[3][3]-----------\n");
	printf(" (0:gas, 1:fluid, 2: solid )\n");
	for(i=0; i<3; i++){
	for(j=0; j<3; j++){
		printf("%lf ",gmm[i][j]);
	}
	printf("\n");
	}
};


Cell::Cell(){
	for(int i=0;i<8;i++){
		cnct[i]=0;
	};
	ID=-1;
	iad=-1;
	erg=0.0;
	erg_bff=0.0;

};
PoreCells::PoreCells(){
	ncell=0;
};
void PoreCells::setup(
	Ellip *els,int nelp,	// ellipses
	bool set_opr, 		// set operator true/false=union/intersection
	int Lev_Max,		// level
	Bbox bx			// bounding box
){
	Tree4::setup(els,nelp,set_opr,Lev_Max,bx); // generate quad tree to manage the domain
	Tree4::set_grid_params(); // set regular grid parameters
	Tree4::draw();

	double *xc,xf[2],xg[2];
	int i,j,k,l,m;
	int ityp,jtyp,isum,iad;

	ncell=0;
	for(i=0;i<n_leaves;i++){
	       	if(leaves[i].isin()>0) ncell++;	// count number of pore cells
	};
	printf("ncell=%d\n",ncell);
	cells=(Cell *)malloc(sizeof(Cell)*ncell); // allocate memory
	for(i=0;i<ncell;i++){
		cells[i].phs=0;	// all cells set to gas phase
		cells[i].phs_bff=0;
	};

	isum=0;
	iad=0;
	int ix,iy;
	for(i=0;i<Nx;i++){
		xf[0]=Xa[0]+dx[0]*(i+0.5);
	for(j=0;j<Ny;j++){
		xf[1]=Xa[1]+dx[1]*(j+0.5);
		ityp=QtreeFind(&qp0,xf);
		if(ityp>0){
			cells[iad].ID=isum;	// linear grid index 
			cells[iad].iad=iad;	// data address in cells[iad];
			cells[iad].bnd=false;
			if(ityp==1) cells[iad].bnd=true;
			iad++;
		}	// end if
		isum++;
	}	// end_j
	}	// end_i
	printf("iad=%d\n",iad);
};
void PoreCells::l2ij(int l, int *i, int *j){
	(*i)=l/Ny;
	(*j)=l%Ny;
};
int PoreCells::find(int id){

	if(id < cells[0].ID) return(-1);
	if(id > cells[ncell-1].ID) return(-1);
	int i1=0;
	int i2=ncell-1;
	int im;

	if(cells[i1].ID==id) return(i1);
	if(cells[i2].ID==id) return(i2);

	while(i2-i1>1){
		im=floor((i1+i2)*0.5);
		if(cells[im].ID==id) return(im);
		if(cells[im].ID>id){
			i2=im;
		}else{
			i1=im;
		}
	};
	if(cells[i1].ID==id) return(i1);
	if(cells[i2].ID==id) return(i2);
	return(-1);
};

void PoreCells::connect(){
	FILE *fp=fopen("temp.out","w");
	int ic;
	int i,j,k,l,m;
	int ix,iy,id,iad;
	for(ic=0; ic<ncell;ic++){
		l2ij(cells[ic].ID,&i,&j);
		m=0;
		cells[ic].nc=0;
		for(k=-1; k<=1; k++){
			ix=i+k;
			if(ix<0) ix+=Nx;
			if(ix>=Nx) ix-=Nx;
		for(l=-1; l<=1; l++){
			if((k*k+l*l)==0) continue;
			iy=j+l;
			if(iy<0) iy+=Ny;
			if(iy>=Ny) iy-=Ny;
			id=ix*Ny+iy;
			iad=find(id);
			if(iad!=-1){
				cells[ic].cnct[m]=iad;
				cells[ic].nc++;
				cells[ic].cncl[m]=cells+iad;
				m++;
			};
		}
		}
		fprintf(fp,"cell no.=%d, ID=%d\n",ic,cells[ic].ID);
		fprintf(fp," number of connected nodes=%d\n",cells[ic].nc);
	};
	fclose(fp);
};

int PoreCells::init(double sr){

	std::mt19937_64 engine(-2);
	std::uniform_int_distribution<int>MT01(0,ncell);

	Sr=sr;
	n_water=int(ncell*Sr);
	n_void=ncell-n_water;
	indx_w=(int *)malloc(sizeof(int)*n_water);
	indx_v=(int *)malloc(sizeof(int)*n_void);

	int count=0,ii;
	while(count<=n_water){
		ii=MT01(engine);
		if(cells[ii].phs==0){
		       cells[ii].phs=1;	// set to fluid phase
		       cells[ii].phs_bff=1;
		       count++;
		};
	}

	int iw=0,iv=0;
	for(int i=0;i<ncell;i++){
		if(cells[i].phs==0) indx_v[iv++]=i;
		if(cells[i].phs==1) indx_w[iw++]=i;
	};
	printf("ncell=%d\n",ncell);
	printf("n_water=%d, n_void=%d\n",n_water,n_void);
	printf("iw=%d, iv=%d\n",iw,iv);
	return(n_water);
};
void PoreCells::load_gmm(
	double th	// contact angle in degree
){
	mtrl.load(th);
	mtrl.print_gmm();
	mtrl.normalize();
	mtrl.print_gmm();
};
double PoreCells::cell_energy(int iad){

	Cell *cell_i,*cell_j;

	int iphs,jphs;
	int k,nc;
	double Erg;

	cell_i=cells+iad;
	nc=cell_i->nc;
	iphs=cell_i->phs;
	Erg=0.0;
	for(k=0;k<nc;k++){
		cell_j=cell_i->cncl[k];
		jphs=cell_j->phs;
		Erg+=mtrl.gmm[iphs][jphs];
	}
	Erg+=(8-nc)*mtrl.gmm[iphs][2];
	return(Erg*0.5);
};
double PoreCells::total_energy(){
	Etot=0.0;
	double dE;
	for(int i=0;i<ncell;i++){
		dE=cell_energy(i);
	       	Etot+=dE;
		cells[i].erg=dE;
	}
	return(Etot);
};
double PoreCells::swap(int id, int jd){

	cells[id].phs_bff=cells[id].phs;
	cells[id].erg_bff=cells[id].erg;

	cells[jd].phs_bff=cells[jd].phs;
	cells[jd].erg_bff=cells[jd].erg;

	int itmp;
	itmp=cells[id].phs;
	cells[id].phs=cells[jd].phs;
	cells[jd].phs=itmp;

	double Ei,Ej,dEi,dEj;
	Ei=cell_energy(id);
	Ej=cell_energy(jd);

	dEi=Ei-cells[id].erg;
	dEj=Ej-cells[jd].erg;

	cells[id].erg=Ei;
	cells[jd].erg=Ej;

	return(dEi+dEj);
	
};
void PoreCells::reject_swap(int id, int jd){

	cells[id].phs=cells[id].phs_bff;
	cells[id].erg=cells[id].erg_bff;
	cells[jd].phs=cells[jd].phs_bff;
	cells[jd].erg=cells[jd].erg_bff;

};

double PoreCells::MC_stepping(Temp_Hist TH){
	int i,j,itmp;
	int id,jd;
	double dE,prb;
	static std::mt19937_64 engine(-2);
	std::uniform_int_distribution<int>MTv(0,n_void-1);
	std::uniform_int_distribution<int>MTw(0,n_water-1);
	std::uniform_real_distribution<double>Urnd(0.0,1.0);
	double E0=Etot;

	for(i=0;i<n_water;i++){
		id=indx_w[i];
		j=MTv(engine);
		jd=indx_v[j];
		dE=swap(id,jd);
		prb=exp(-dE/TH.Temp);
		if(Urnd(engine)<=prb){
			itmp=indx_w[i];
			indx_w[i]=indx_v[j];
			indx_v[j]=itmp;
			Etot+=dE;
		}else{
			reject_swap(id,jd);
		};
	};

	for(i=0;i<n_void;i++){
		id=indx_v[i];
		j=MTw(engine);
		jd=indx_w[j];
		dE=swap(id,jd);
		prb=exp(-dE/TH.Temp);
		if(Urnd(engine)<=prb){
			itmp=indx_v[i];
			indx_v[i]=indx_w[j];
			indx_w[j]=itmp;
			Etot+=dE;
		}else{
			reject_swap(id,jd);
		};
	};
	return(Etot-E0);
};
void PoreCells::write_phs(){

	FILE *fp;
	int i,ix,iy,ID;
	double xx,yy;

	fp=fopen("kcell_w.out","w");
	for(i=0;i<n_water;i++){
		ID=cells[indx_w[i]].ID;
		l2ij(ID,&ix,&iy);
		xx=Xa[0]+(ix+0.5)*dx[0];
		yy=Xa[1]+(iy+0.5)*dx[1];
		fprintf(fp,"%lf %lf\n",xx,yy);
	};
	fclose(fp);

	fp=fopen("kcell_v.out","w");
	for(i=0;i<n_void;i++){
		ID=cells[indx_v[i]].ID;
		l2ij(ID,&ix,&iy);
		xx=Xa[0]+(ix+0.5)*dx[0];
		yy=Xa[1]+(iy+0.5)*dx[1];
		fprintf(fp,"%lf %lf\n",xx,yy);
	};
	fclose(fp);
};
void PoreCells::fwrite_cells(char fn[128]){
	double PI=4.0*atan(1.0);
	FILE *fp=fopen(fn,"w");
	fprintf(fp,"# thE (contact angle [deg])\n");
	fprintf(fp,"%lf\n",mtrl.thE/PI*180.0);
	fprintf(fp,"# Sr (degree of saturation)\n");
	fprintf(fp,"%lf\n",Sr);

	fprintf(fp,"# Xa[2], Xb[2]\n");
	fprintf(fp,"%lf, %lf\n",Xa[0],Xa[1]);
	fprintf(fp,"%lf, %lf\n",Xb[0],Xb[1]);
	fprintf(fp,"# Nx, Ny\n");
	fprintf(fp,"%d, %d\n",Nx,Ny);
	fprintf(fp,"# ncell (number of fluid + gas cells)\n");
	fprintf(fp,"%d\n",ncell);
	fprintf(fp,"# Global cell No. (ID), phase(0: gas, 1:fluid), boundary cell (1:True,0:False)\n");
	for(int iad=0;iad<ncell;iad++){
		fprintf(fp,"%d %d %d\n",cells[iad].ID, cells[iad].phs,cells[iad].bnd);
	};
	fclose(fp);
};
int PoreCells::grid_type(int i, int j){

	int k,l,ix,iy,nc=0;
	int iad,phs,ityp;
	int ngrid[3]={0,0,0}; // solid, fluid, gas

	for(k=-1; k<1; k++){	
		ix=i+k;
		if(ix<0) ix+=Nx;
		if(ix>=Nx) ix-=Nx;
	for(l=-1; l<1;l++){	
		iy=j+l;
		if(iy<0) iy+=Ny;
		if(iy>=Ny) iy-=Ny;
		l=ix*Ny+iy;
		iad=find(l);

		if(iad==-1){
			ngrid[2]++;
		}else{
			phs=cells[iad].phs;
			ngrid[phs]++;
		}
	}
	}

	ityp=2;	// solid 
	if(ngrid[1]>0){
	     ityp=1;	// fluid
	}else if(ngrid[0]>0){
	     ityp=0;	// gas
	}

	return(ityp);	// 0:gas, 1:fluid, 2:solid
};
int PoreCells::count_grids(){
	int i,j,ityp,jtyp;
	int ng=0;
	for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
		ityp=PoreCells::grid_type(i,j);
		jtyp=PoreCells::grid_loc(i,j);
		if((ityp==1) || (jtyp==1)) ng++;
//		if( ityp==1 ) ng++;
	}
	}
	return(ng);
};
int PoreCells::grid_loc(int i, int j){

	int k,l,ix,iy;
	int iad,ityp;
	int ngrid[3]={0,0,0}; // int, bnd, ext (relative to solid phase cell)

	for(k=-1; k<1; k++){	
		ix=i+k;
		if(ix<0) ix+=Nx;
		if(ix>=Nx) ix-=Nx;
	for(l=-1; l<1;l++){	
		iy=j+l;
		if(iy<0) iy+=Ny;
		if(iy>=Ny) iy-=Ny;
		l=ix*Ny+iy;
		iad=find(l);

		if(iad==-1){
			ngrid[0]++;	// interior point (solid)
		}else if(cells[iad].bnd){
			ngrid[1]++;	// boundary point (solid/pore)
		}else{
			ngrid[2]++;	// exterior point (pore)
		}
	}
	}

	ityp=0;	// interior point (solid)
	if(ngrid[1]>0){
	     ityp=1;	// boundary point (solid/pore interface) 
	}else if(ngrid[0]>0){
	     ityp=0;	// exterior point (pore) 
	}

	return(ityp);	// 0:interior, 1:interfacial , 2:exterior points 
};
void PoreCells::load_cell_data(char fn[128]){
	FILE *fp=fopen(fn,"r");
	char cbff[128];
	double th;

	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&th);
	load_gmm(th);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&Sr);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf, %lf\n",Xa,Xa+1);
	fscanf(fp,"%lf, %lf\n",Xb,Xb+1);

	fgets(cbff,128,fp);
	fscanf(fp,"%d, %d\n",&Nx,&Ny);

	Wd[0]=Xb[0]-Xa[0];
	Wd[1]=Xb[1]-Xa[1];
	dx[0]=Wd[0]/Nx;
	dx[1]=Wd[1]/Ny;

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&ncell);
	fgets(cbff,128,fp);

	cells=(Cell *)malloc(sizeof(Cell)*ncell); // allocate memory

	int ID,phs,bnd;
	for(int i=0;i<ncell;i++){
		fscanf(fp,"%d %d %d\n",&ID, &phs, &bnd);
		cells[i].ID=ID;
		cells[i].phs=phs;
		cells[i].bnd=true;
		if(bnd==0) cells[i].bnd=false;
	};


	fclose(fp);
};
void PoreCells::write_leaves(){
	Tree4::draw();
};
