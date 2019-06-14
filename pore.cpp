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
/*
//---------------------------------------------------------------------------------
class Material{
	public:
		double gmm[3][3];	// interfacial energy (0:gas, 1:fluid, 2: solid)
		double thE;		// contact angle
		void load(double th); 	// set gmm values 
		double gmm_max;		// maximum value in gmm
		void  normalize();	// normalize gmm
		void  print_gmm(); // print interfacial energy
		double erg;	// cell energy
};
class Cell{
	public:
		QPatch *qp0;
		int cnct[8];//connected ?
		Cell *cncl[8]; //pointer to connected cells
		int nc;	// number  of connected cells
		Cell();	// constructuor
		int ID;		// linear index for 2D Grid
		int iad;	// address in PoreCells[ncell]
		int phs;	// 0=gas, 1=fluid, 2=solid 
		int phs_bff;
		double erg;
		double erg_bff;
	private:
};
class PoreCells:public Tree4{
	public:
		int ncell; 	// number of cells
		Cell *cells;	// pointer(array) to cell class instances
		PoreCells();	// constructor
		void setup(Ellip *els,int nelp,bool set_opr, int Lev_Max,Bbox bx); // generate cells 
		int find(int id);	// find cell having a given linear grid number(id).
		void connect();	// establish neghboring cell connection 
		void l2ij(int l, int *i, int *j); // index transform ( linear to 2D index)
		double Sr;	// degree of water saturation
		int n_void,n_water; // number of gas and fluid cells, resp.
		int init(double sr); // initialize phase distribution
		int *indx_w;	// index set of fluid phase
		int *indx_v;	// index set of gas phase
		double cell_energy(int iad); // evaluate interfacial energy/cell
		double total_energy(); // evaluate total interfacial energy
		Material mtrl;	// material constants (interfaceial energy)
		void load_gmm(double th); // load gmm data (th = fluid/solid contact angle)
		double Etot;	// total interfacial energy
		double swap(int id, int jd); // swap cell id & jd tempralily
		void reject_swap(int id, int jd); // apply swap 
		double MC_stepping(Temp_Hist TH);	// Monte Carlo stepping 
		void write_phs();
		int count_grids();
		int grid_type(int i,int j);
	private:
};
*/
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
			/*-------------------------------------------
			cells[iad].nc=0;
			k=0;
			for(l=-1;l<=1;l++){
				xg[0]=xf[0]+l*dx[0];
				ix=i+l;
				if(ix<0) ix+=Nx;
				if(ix>=Nx) ix-=Nx;
			for(m=-1;m<=1;m++){
				if((l*l+m*m)==0) continue;
				iy=j+m;
				if(iy<0) iy+=Ny;
				if(iy>=Ny) iy-=Ny;

				xg[1]=xf[1]+m*dx[1];
				jtyp=QtreeFind(&qp0,xg);
				if(jtyp>0){
					cells[iad].cnct[k]=ix*Ny+iy;	// linear grid index
					cells[iad].nc++;
					k++;
				}
			}	// end_m
			}	// end_l
			*/
			//-------------------------------------------
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

void PoreCells::connect0(){
	FILE *fp=fopen("temp0.out","w");
	int i,j,iad;
	int ix,iy;
	int jx,jy;
	int nrm,ixd,iyd;
	for(i=0;i<ncell;i++){
		fprintf(fp,"cell no.=%d, ID=%d\n",i,cells[i].ID);
		fprintf(fp," number of connected nodes=%d\n",cells[i].nc);
		l2ij(cells[i].ID,&ix,&iy);
		for(j=0;j<cells[i].nc;j++){
			iad=find(cells[i].cnct[j]);
			cells[i].cncl[j]=cells+iad;
			l2ij(cells[i].cncl[j]->ID,&jx,&jy);
			ixd=abs(ix-jx);
			if(Nx-ixd < ixd) ixd=Nx-ixd;
			iyd=abs(iy-jy);
			if(Ny-iyd < iyd) iyd=Ny-iyd;
			nrm=ixd+iyd;
			fprintf(fp," %d",nrm);
			if(nrm>2) printf("connectivity error !");
		};
		fprintf(fp,"\n");
	}
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
	std::uniform_int_distribution<int>MTv(0,n_void);
	std::uniform_int_distribution<int>MTw(0,n_water);
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

	return(ityp);
};
int PoreCells::count_grids(){
	int i,j,ityp;
	int ng=0;
	for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
		ityp=PoreCells::grid_type(i,j);
		if(ityp==1) ng++;
	}
	}
	return(ng);
};
void PoreCells::load_cell_data(char fn[128]){
	FILE *fp=fopen(fn,"r");
	char cbff[128];

	double th;
//	double Xa[2],Xb[2];

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
/*
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
	//Pcll.connect(); // establish connection among pore coverning cells
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
	return(0);
};
*/
