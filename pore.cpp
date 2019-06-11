#include<stdio.h>
#include<stdlib.h>
#include<random>
#include"set2d.h"
#include"vec2.h"
#ifndef __TCNTRL__
	#define __TCNTRL__
	#include "tcntrl.h"
#endif

using namespace std;
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
		int phs;	// 0 = solid, 1 = fluid, 2 = void
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
	private:
};
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
		cells[i].phs=0;
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
			cells[iad].nc=0;
			//-------------------------------------------
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

	while(i2-i1>0){
		im=(i1+i2)*0.5;
		if(cells[im].ID==id) return(im);
		if(cells[im].ID>id){
			i2=im;
		}else{
			i1=im;
		}
	};
	return(-1);
};

void PoreCells::connect(){
	FILE *fp=fopen("temp.out","w");
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
			//fprintf(fp," %d (%d)",cells[i].cnct[j],iad);
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
		       cells[ii].phs=1;
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

	//printf("Etot(initial)=%lf\n",E0);
	//printf("Etot(final  )=%lf\n",Etot);
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
	printf("Etot=%lf\n",Etot);
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

int main(){


	Solid sld;

	int Lev=9;	// Quad-tree height
	double Wd[2]={1.0,1.0}; // Unit Cell Size
	double Xa[2]={0.0,0.0}; // Unit Cell position (lowerleft vertex)
	char fn[128]="solid.dat";	// solid phase data file

	puts(fn);
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
	int nstep=300;
	Temp_Hist TH(T1,T2,nstep);

	double dE;
	while(TH.cont_iteration){
		dE=Pcll.MC_stepping(TH);
		printf("T=%lf, dE=%lf\n",TH.Temp,dE);
		TH.inc_Temp_exp();
	};
	Pcll.write_phs();
	//Pcll.draw();
	return(0);
};
