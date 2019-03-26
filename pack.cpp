#define DB 0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>
#include "set2d.h"
#ifndef __TCNTRL__
	#define _TCNTRL
	#include "tcntrl.h"
#endif

using namespace std;

class Node{
	public:
		int iad;// data count (address in the Node class array NDs[]) 
		int id;	// grid number in the underlying regular grid
		int cnct[4];	// connected  nodes
		Node *cnds[4];	// connected nodes 
		Node();	//default constructor
		int nc;
	private:
};
Node::Node(){
	id=-1;
	for(int i=0;i<4;i++){
	       	cnct[i]=-1;
		cnds[i]=NULL;
	}
};

class Grid{
	public:
		int Ng;		// total number of grids
		int Nx,Ny;	// number of underling regular grid (Ng != Nx x Ny)
		Node *NDs;	// grid points
		bool ready;
		double Xa[2],Wd[2],Xb[2],dx[2];
		int Ndiv[2];
		Grid();	// default constructor
		void setup(Tree4 tr4);
		int find(int  k);
		void l2cod(int l,double *x, double *y);
		void indx2cod(int i,int j, double *x, double *y);
		void l2ij(int l, int *i, int *j);
		void grid_cod(int inod,double *xcod, double *ycod);
	private:
		void mem_alloc();
};

int Grid::find(int k){
	if(k<NDs[0].id) return(-1);
	if(k>NDs[Ng-1].id) return(-1);
	int i1=0;
	int i2=Ng-1;
	int im;

	if(NDs[i1].id==k) return(i1);
	if(NDs[i2].id==k) return(i2);
	while(i2-i1>0){
		im=(i1+i2)/2;
		if(NDs[im].id==k) return(im);
		if(NDs[im].id>k){
			i2=im;
		}else{
			i1=im;
		}
	}
	return(-1);
};

void Grid::setup(Tree4 tr4){
	tr4.set_grid_params();

	Nx=tr4.Nx;
	Ny=tr4.Ny;
	for(int i=0;i<2;i++){
		Xa[i]=tr4.Xa[i];
		Xb[i]=tr4.Xb[i];
		Wd[i]=tr4.Wd[i];
		dx[i]=tr4.dx[i];
	}

	int cnct[4];
	int ngb=0;
	int m,k;
	for(k=0; k<tr4.Nx; k++){
	for(m=0; m<tr4.Ny; m++){
		if(tr4.grid_type(k,m,cnct)==1) ngb++;
	}
	}
	printf("ngb=%d\n",ngb);
	Ng=ngb;
	Grid::mem_alloc();

	ngb=0;
	int isum=0,l,nc;
	for(k=0; k<tr4.Nx; k++){
	for(m=0; m<tr4.Ny; m++){
		if(tr4.grid_type(k,m,cnct)==1){
			NDs[ngb].id=isum;
			NDs[ngb].iad=ngb;
			nc=0;
			for(l=0;l<4;l++){
			       NDs[ngb].cnct[l]=cnct[l];
			       nc+=cnct[l];
			};
			NDs[ngb].nc=nc;
			ngb++;
		};
		isum++;
	}
	}

	int iofst[4]={-1, 0, 1, 0};
	int jofst[4]={ 0,-1, 0, 1};
	int ix0,jy0,ix,jy,adr,id,ig;
	for(k=0;k<Ng;k++){
		ix0=NDs[k].id/Ny;
		jy0=NDs[k].id%Ny;
		nc=0;
		for(l=0;l<4;l++){
			if(NDs[k].cnct[l]==1){
				ix=ix0+iofst[l];		
				jy=jy0+jofst[l];		
				while(ix<0) ix+=Nx;
				while(jy<0) jy+=Ny;
				while(ix>=Nx) ix-=Nx;
				while(jy>=Ny) jy-=Ny;
				id=ix*Ny+jy;
				adr=Grid::find(id);
				NDs[k].cnds[nc]=(NDs+adr);
				nc++;
			};
		};
	};

};
void Grid::l2cod(int l,double *x,double *y){
	int i,j;
	Grid::l2ij(l,&i,&j);
	Grid::indx2cod(i,j,x,y);
};
void Grid::indx2cod(int i,int j,double *x,double *y){
	*x=Xa[0]+dx[0]*i;
	*y=Xa[1]+dx[1]*j;
};
void Grid::l2ij(int l, int *i, int *j){
	*i=l/Ny;
	*j=l%Ny;
};
Grid::Grid(){
	Ng=1;
	ready=false;
};
void Grid::mem_alloc(){
	NDs=(Node *)malloc(sizeof(Node)*Ng);
	ready=true;
};
void Grid::grid_cod(int inod,double *xcod, double *ycod){
	int id=NDs[inod].id;
	int i,j;
	Grid::l2ij(id,&i,&j);
	Grid::indx2cod(i,j,xcod,ycod);
};
//----------------------------------------------------------
#if DB==4
// Testing Grid and Node classes on which random walks will be taken
int main(){
	int np=100;	// number of particles
	int Lev=9;	// Quad-tree height
	double Wd[2]={1.0,1.0}; // Unit Cell Size
	//Solid sld(np,Wd);
	char fn[128]="solid.dat";
	Solid sld;
	
	double Xa[2]={0.0,0.0};
	sld.load(fn);
	sld.bbox.setup(Xa,Wd);

	Tree4 tr4;
	tr4.setup(sld.els,sld.nelp,false,Lev,sld.bbox);

	double xf[2]={0.5,0.384};
	printf("Is interior point? --> %d\n",QtreeFind(&(tr4.qp0),xf));
	tr4.count();
	tr4.set_grid_params();
	printf("Xa=%lf %lf\n",tr4.Xa[0],tr4.Xa[1]);
	printf("Xb=%lf %lf\n",tr4.Xb[0],tr4.Xb[1]);
	printf("Wd=%lf %lf\n",tr4.Wd[0],tr4.Wd[1]);
	printf("dx=%lf %lf\n",tr4.dx[0],tr4.dx[1]);

	Grid gd;
	gd.setup(tr4);


	std::mt19937_64 engine(-2);
	std::uniform_real_distribution<double>MT01(0.0,1.0);
	double xcod,ycod;
	int next,now;
	printf("Start grid=%d\n",now);
	
	int nwk=10000,Nt=20000,inc=100; 
	
	//Node *nd0[200];
	Node **nd0=(Node **)malloc(sizeof(Node*)*nwk);
	int i,j;
	FILE *fp=fopen("rw.out","w");
	for(i=0;i<nwk;i++){
		now=int(MT01(engine)*gd.Ng);
		nd0[i]=&gd.NDs[now];
	};

	fprintf(fp,"%d,%d,%d\n",nwk,Nt,inc);
	fprintf(fp,"%le,%le\n",tr4.dx[0],tr4.dx[1]);
	double x0,y0;
	double tolx=tr4.dx[0]*1.001;
	double toly=tr4.dx[1]*1.001;
	double *ofx=(double *)malloc(sizeof(double)*nwk);
	double *ofy=(double *)malloc(sizeof(double)*nwk);
	printf("tols=%lf %lf\n",tolx,toly);
	for(j=0;j<Nt;j++){
		for(i=0;i<nwk;i++){
			gd.grid_cod(nd0[i]->iad,&x0,&y0);

			next=int(MT01(engine)*nd0[i]->nc);
			nd0[i]=nd0[i]->cnds[next];
			gd.grid_cod(nd0[i]->iad,&xcod,&ycod);

			if(xcod-x0>tolx) ofx[i]-=Wd[0];
			if(x0-xcod>tolx) ofx[i]+=Wd[0];
			if(ycod-y0>toly) ofy[i]-=Wd[1];
			if(y0-ycod>toly) ofy[i]+=Wd[1];
			if(j%inc==0) fprintf(fp,"%lf,%lf\n",xcod+ofx[i],ycod+ofy[i]);
		}
	}


	/*
	int i,ig;
	for(i=0;i<10;i++){
		ig=gd.NDs[i].id;
		printf("grid no.=%d, address=%d\n",ig,gd.find(ig));
	};
	*/
	tr4.draw();
	return(0);
};
#endif

#if DB==3	// Debugging memory leak 
// -->  memory blocks for Poly::xs,ys are not freed 
int main(){
	int np=200;
	double Wd[2]={1.0,1.0};
	char fn[128]="geom0.dat";

	Solid sld(np,Wd);
	sld.draw(fn,50);

	Tree4 tr;
	for(int i=0;i<1;i++){
		printf("i_try=%d\n",i);
		tr.setup(sld.els,sld.nelp,false,10,sld.bbox);
//		tr.setup(sld,10);
		tr.draw();
		tr.clean();
	};

	int lev_max=9;
	printf("S=%lf\n",sld.area(lev_max));


/*
	Tree4 tr4;
	tr4.setup(sld,9);
	tr4.draw();
	tr4.clean();
*/
};
#endif

#if DB==0
int main(){
	//std::uniform_real_distribution<double> RndR;
	//RndR=uniform_real_distribution<double>(0.0,1.0);
	Temp_Hist TH;
	char fnt[128]="temp_hist.inp";	// Annealing parameter (temperature control)
	TH.load(fnt);

	int np=400;	// number of particles
	int Lev=9;	// Quad-tree height
	double Wd[2]={1.0,1.0}; // Unit Cell Size
	char fn[128]="geom0.dat";
	double dE_tot=0.0,dE;

	Solid sld(np,Wd);
	sld.write(fn);
	puts("initial geometry written\n");
	int tmp;
	scanf("%d",&tmp);
	//sld.draw(fn,50);

	double Xa[2]={0.0,0.0};
	char fns[128]="solid1.dat";
	/*
	Solid sld;
	sld.load(fns);
	sld.bbox.setup(Xa,Wd);
	*/

	sld.area(Lev);
	FILE *fl=fopen("log.dat","w");
	while(TH.cont_iteration){
		dE=sld.MC(TH);
		dE_tot+=dE;
		TH.inc_Temp_exp();
		printf("tau=%lf, dE=%le dE=%le\n",TH.tau(),dE,dE_tot);
		fprintf(fl,"%ld %le %le %le %le\n",TH.istep,TH.tau(),TH.Temp,dE,dE_tot);
		fflush(stdout);
	};
	printf("dE_tot=%lf\n",dE_tot);
	sprintf(fn,"geom1.dat");
	sld.draw(fn,50);
	sld.area(Lev);

	sprintf(fn,"solid.dat");
	sld.write(fn);

	Tree4 tr4;
	tr4.setup(sld.els,sld.nelp,false,Lev,sld.bbox);
	tr4.draw();
	tr4.clean();
	return(0);
};
#endif


#if DB==1	
// checking area(el1,el2,levmax,isect)
// Introducing Tree4 class with which quad-tree based quadrature can be done 
int main(){

	Ellip el1,el2;
	el1.set_xc(0.0,-0.0);
	el1.set_radi(1.0,0.5);
	el1.set_phi(0.0);
	el1.set_bbox();
	printf("A1=%lf\n",el1.area());

	el2.set_xc(0.0,0.5);
	el2.set_radi(1.0,0.5);
	el2.set_phi(0.0);
	el2.set_bbox();
	printf("A2=%lf\n",el2.area());

	bool isect=true;
	int lev_max=6;
	double S=area(el1,el2,lev_max,isect);
	printf("A=%lf, err=%lf\n",S,S-el1.area()-el2.area());

	Tree4 tr4;
	tr4.setup(el1,el2,isect,lev_max);
	tr4.draw();
	printf("A=%lf\n",tr4.area());
	tr4.clean();

	char fname[128]="geom.dat",md[3]="w";	
	el1.draw(fname,100,md);
	sprintf(md,"a");	
	el2.draw(fname,100,md);


	return(0);
};
#endif
#if DB==2
int main(){
	std::random_device rd;
	//std::mt19937 mt(rd());
	std::mt19937_64 mt(-1);
	std::uniform_real_distribution<double> RndR(0.0,1.0);

	char fn[128]="log.txt";
	char md[3]="w",mode[3];

	double PI=4.0*atan(1.0);
	int i,np=100;
	double x,y;
	double ra,rb,aspect=0.6;
	double phi;
	double Wd[2];

	Solid sld(np);
	Wd[0]=1.0;
	Wd[1]=1.0;

	for(i=0;i<np;i++){
		x=RndR(mt)*Wd[0];	
		y=RndR(mt)*Wd[1];	
		ra=0.1*Wd[0];
		rb=ra*aspect;
		phi=RndR(mt)*PI;
		//printf("x,y=%lf %lf, phi=%lf[deg]\n",x,y,phi/PI*180.0);
		sld.els[i].set_xc(x,y);
		sld.els[i].phi=phi;
		sld.els[i].set_radi(ra,rb);
		sld.isect[i]=false;
		sld.els[i].set_bbox();
	};

	sprintf(fn,"geom.dat");
	sld.draw(fn,50);

	int q,p=int(RndR(mt)*np);
	bool pq;
	Ellip elp=sld.els[p],elq;
	//elp.draw(50);
	//ep.bbox.draw();
	double dE=0.0,S=0.0;
	Tree4 tr4;
	for(q=0;q<np;q++){
		if(q==p) continue;
		pq=bbox_cross(elp,sld.els[q]);
		if(pq){
			elq=sld.els[q];
		       	//elq.draw(50);
		       	//sld.els[q].bbox.draw();
			tr4.setup(elp,elq,true,6);
			S=tr4.area();
			tr4.clean();
			printf("S=%lf\n",S);
			dE+=S;
		}
	}	
	printf("dE=%lf\n",dE);

	double Xa[2],Xb[2];
	Xa[0]=0.0; 
	Xa[1]=0.0;
	Xb[0]=Xa[0]+Wd[0];
	Xb[1]=Xa[1]+Wd[1];
	sld.bbox.setup(Xa,Xb);

	Ellip elr=elp;
	printf("xc=%lf %lf -->",elr.xc[0],elr.xc[1]);
	double ux,uy;
	ux=RndR(mt)-0.45;
	uy=RndR(mt)-0.15;
	elr.slide(ux,uy,sld.bbox);
	printf("xc=%lf %lf\n",elr.xc[0],elr.xc[1]);

	dE=0.0; S=0.0;
	for(q=0;q<np;q++){
		if(q==p) continue;
		pq=bbox_cross(elr,sld.els[q]);
		if(pq){
			elq=sld.els[q];
			tr4.setup(elr,elq,true,6);
			S=tr4.area();
			tr4.clean();
			printf("S=%lf\n",S);
			dE+=S;
		}
	}	
	printf("dE=%lf\n",dE);

	sld.perturb(p,ux,uy,0.0);

	return(0);
};
#endif
