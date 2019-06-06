#include<stdio.h>
#include<stdlib.h>
#include<random>
#include"set2d.h"
#include"vec2.h"

using namespace std;
//---------------------------------------------------------------------------------
class Cell{
	public:
		QPatch *qp0;
		int cnct[8];//connected ?
		Cell *cncl[8];
		int nc;	// number 
		Cell();
		int ID;		// linear index for 2D Grid
		int iad;	// address in PoreCells[ncell]
		int phs;	// 0 = solid, 1 = fluid, 2 = void
	private:
};
class PoreCells:public Tree4{
	public:
		int ncell; 
		Cell *cells;
		PoreCells();
		void setup(Ellip *els,int nelp,bool set_opr, int Lev_Max,Bbox bx);
		int find(int id);
		void connect();
		void l2ij(int l, int *i, int *j);
	private:
};
//---------------------------------------------------------------------------------
Cell::Cell(){
	for(int i=0;i<8;i++){
		cnct[i]=0;
	};
	ID=-1;
	iad=-1;
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
int main(){

	printf("start program\n");
	Solid sld;

	int Lev=9;	// Quad-tree height
	double Wd[2]={1.0,1.0}; // Unit Cell Size
	double Xa[2]={0.0,0.0}; // Unit Cell position (lowerleft vertex)
	char fn[128]="solid.dat";	// solid phase data file

	puts(fn);
	sld.load(fn);	// import solid phase data
	sld.bbox.setup(Xa,Wd); // set bounding box

	PoreCells Pcll;
	Pcll.qp0.refine[0]=true;
	Pcll.qp0.show_prms();
	Pcll.setup(sld.els,sld.nelp,false,Lev,sld.bbox);
	Pcll.connect();

	//Pcll.draw();
	return(0);
};
