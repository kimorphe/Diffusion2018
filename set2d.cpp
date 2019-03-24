#define DB 0
#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
#include "vec2.h"
#include "set2d.h"

using namespace std;
void Bbox::set_Xa(double x, double y){
	Xa[0]=x; Xa[1]=y;
};
void Bbox::set_Xb(double x, double y){
	Xb[0]=x; Xb[1]=y;
};
void Bbox::set_Wd(){
	Wd[0]=Xb[0]-Xa[0];
	Wd[1]=Xb[1]-Xa[1];
};
void Bbox::setup(double xa[2], double xb[2]){
	for(int i=0;i<2;i++){
		Xa[i]=xa[i];
		Xb[i]=xb[i];
		Wd[i]=Xb[i]-Xa[i];
	}
};
double Bbox::area(){
	return(Wd[0]*Wd[1]);
};
void Bbox::draw(char fn[128],char mode[3]){
	FILE *fp=fopen(fn,mode);

	fprintf(fp,"%lf %lf\n",Xa[0],Xa[1]);
	fprintf(fp,"%lf %lf\n",Xb[0],Xa[1]);
	fprintf(fp,"%lf %lf\n",Xb[0],Xb[1]);
	fprintf(fp,"%lf %lf\n",Xa[0],Xb[1]);
	fprintf(fp,"%lf %lf\n",Xa[0],Xa[1]);

	fprintf(fp,"\n");
	fclose(fp);
};
void Bbox::draw(){

	printf("%lf %lf\n",Xa[0],Xa[1]);
	printf("%lf %lf\n",Xb[0],Xa[1]);
	printf("%lf %lf\n",Xb[0],Xb[1]);
	printf("%lf %lf\n",Xa[0],Xb[1]);
	printf("%lf %lf\n",Xa[0],Xa[1]);
	puts("");
};
void Bbox::slide(double ux, double uy){
	Xa[0]+=ux; Xb[0]+=ux; 
	Xa[1]+=uy; Xb[1]+=uy; 
};
double dmin(double x, double y){
	if(x<=y) return(x);
	return(y);
};
double dmax(double x, double y){
	if(x>=y) return(x);
	return(y);
};
Bbox bbox_cross(Bbox b1, Bbox b2){
	Bbox b3;

	double xmin=dmax(b1.Xa[0],b2.Xa[0]);
	double ymin=dmax(b1.Xa[1],b2.Xa[1]);
	double xmax=dmin(b1.Xb[0],b2.Xb[0]);
	double ymax=dmin(b1.Xb[1],b2.Xb[1]);

	if(xmax<xmin) xmax=xmin;
	if(ymax<ymin) ymax=ymin;
	b3.set_Xa(xmin,ymin);
	b3.set_Xb(xmax,ymax);
	b3.set_Wd();
	return(b3);
};
Bbox bbox_union(Bbox b1, Bbox b2){
	Bbox b3;

	double xmin=dmin(b1.Xa[0],b2.Xa[0]);
	double ymin=dmin(b1.Xa[1],b2.Xa[1]);
	double xmax=dmax(b1.Xb[0],b2.Xb[0]);
	double ymax=dmax(b1.Xb[1],b2.Xb[1]);
	b3.set_Xa(xmin,ymin);
	b3.set_Xb(xmax,ymax);
	b3.set_Wd();
	return(b3);
};
bool bbox_cross(Bbox bx, Pixel px){
	double xmin=dmax(bx.Xa[0],px.Xa[0]);
	double ymin=dmax(bx.Xa[1],px.Xa[1]);
	double xmax=dmin(bx.Xb[0],px.Xb[0]);
	double ymax=dmin(bx.Xb[1],px.Xb[1]);
	if(xmax<xmin) return(false);
	if(ymax<ymin) return(false);
	return(true);
};
bool bbox_cross(Ellip el1, Ellip el2){
	double *Xa,*Xb;
	double *Ya,*Yb;
	Xa=el1.bbox.Xa; Xb=el1.bbox.Xb;
	Ya=el2.bbox.Xa; Yb=el2.bbox.Xb;
	double xmin=dmax(Xa[0],Ya[0]);
	double ymin=dmax(Xa[1],Ya[1]);
	double xmax=dmin(Xb[0],Yb[0]);
	double ymax=dmin(Xb[1],Yb[1]);
	if(xmax<xmin) return(false);
	if(ymax<ymin) return(false);
	return(true);
};
//-------------------------------------------------
bool Circ :: isin(double *x){
	double dist;

	dist=(x[0]-xc[0])*(x[0]-xc[0]);
	dist+=(x[1]-xc[1])*(x[1]-xc[1]);
	dist=sqrt(dist);

	if(dist <= radi){
		return true;
	}{
		return false;
	}
};
void Circ::draw(char fn[128],int npnt,char mode[3]){
	FILE *fp=fopen(fn,mode);
	double th1=0.0,th2=atan(1.0)*4.0*2.0;
	double th,dth=(th2-th1)/(npnt-1);
	for(int i=0;i<npnt;i++){
		th=th1+dth*i;
		fprintf(fp,"%lf %lf\n",xc[0]+radi*cos(th),xc[1]+radi*sin(th));
	};
	fprintf(fp,"\n");
	fclose(fp);
};
void Circ::set_bbox(){
	bbox.set_Xa(xc[0]-radi,xc[1]-radi);
	bbox.set_Xb(xc[0]+radi,xc[1]+radi);
};
Pixel::Pixel(){
	Xa[0]=0.0; Xa[1]=0.0;
	Xb[0]=1.0; Xb[1]=1.0;
	ready=false;
};
void Pixel::draw(char fn[128],char mode[3]){
	FILE *fp=fopen(fn,mode);
	if(~ready) Pixel::setup();
	for(int i=0;i<5;i++){
		fprintf(fp,"%lf %lf\n",xs[i%4],ys[i%4]);
	};
	fprintf(fp,"\n");
	fclose(fp);
};
void Pixel::draw(FILE *fp){
	if(~ready) Pixel::setup();
	for(int i=0;i<5;i++){
		fprintf(fp,"%lf %lf\n",xs[i%4],ys[i%4]);
	};
	fprintf(fp,"\n");
};
void Pixel::draw(){
	if(~ready) Pixel::setup();
	for(int i=0;i<5;i++){
		printf("%lf %lf\n",xs[i%4],ys[i%4]);
	};
	printf("\n");
};
void Pixel::print(){
	printf("Xa=%lf %lf\n",Xa[0],Xa[1]);
	printf("Xb=%lf %lf\n",Xb[0],Xb[1]);
	Wd[0]=Xb[0]-Xa[0];
	Wd[1]=Xb[1]-Xa[1];
	Xc[0]=(Xa[0]+Xb[0])*0.5;
	Xc[1]=(Xa[1]+Xb[1])*0.5;
	printf("Wd=%lf %lf\n",Wd[0],Wd[1]);
	printf("Xc=%lf %lf\n",Xc[0],Xc[1]);
};
void Pixel::setup(){
	Wd[0]=Xb[0]-Xa[0];
	Wd[1]=Xb[1]-Xa[1];
	Xc[0]=(Xa[0]+Xb[0])*0.5;
	Xc[1]=(Xa[1]+Xb[1])*0.5;

	xs[0]=Xa[0];
	xs[1]=Xb[0];
	xs[2]=Xb[0];
	xs[3]=Xa[0];

	ys[0]=Xa[1];
	ys[1]=Xa[1];
	ys[2]=Xb[1];
	ys[3]=Xb[1];

	Le[0]=Wd[0];
	Le[1]=Wd[1];
	Le[2]=Wd[0];
	Le[3]=Wd[1];

	ready=true;
};
void Pixel::set_Xa(double xa, double ya){
	Xa[0]=xa;
	Xa[1]=ya;
};
void Pixel::set_Xb(double xb, double yb){
	Xb[0]=xb;
	Xb[1]=yb;
};
double Pixel::dist2pt(double px, double py){
	double x1[2],x2[2];
	int i;
	double r1,r2,r;
	double xf[2],wt[4],hn[4];
	double dx,dy;

	if( ~ready) Pixel::setup();

	wt[0]= px-xs[0];
	wt[2]=-(px-xs[2]);
	wt[1]= py-ys[1];
	wt[3]=-(py-ys[3]);

	hn[0]=-(py-ys[0]);
	hn[2]= (py-ys[2]);
	hn[1]= (px-xs[1]);
	hn[3]=-(px-xs[3]);

	dx=px-xs[0];
	dy=py-ys[0];
	rmin=sqrt(dx*dx+dy*dy);
	rmax=0.0;
	for(i=0;i<4;i++){
		dx=px-xs[i];
		dy=py-ys[i];
		r=sqrt(dx*dx+dy*dy);
		if(r > rmax) rmax=r;
		if(r < rmin) rmin=r;

		if(wt[i]<0.0) continue;
		if(wt[i]>Le[i]) continue;
		if(fabs(hn[i]) < rmin) rmin=fabs(hn[i]);
	};

	is_in=true;
	for(i=0;i<4;i++){
		if(hn[i]>0){
			is_in=false;
			break;
		}
	}
	return(rmin);
};
double Pixel::area(){
	if(!ready) Pixel::setup();
	return(Wd[0]*Wd[1]);	
};
//---------------------------------------------------
int is_cross(Pixel px, Circ cr){
	double rmin,rmax;

	rmin=px.dist2pt(cr.xc[0],cr.xc[1]);
	rmax=px.rmax;
	if(rmax < cr.radi) return(1);	// P < C
	if(rmin > cr.radi){
		if(px.is_in){
			return(2); // C < P
		}else{
			return(0); // C ^ P=\phi
		}
	}
	return(3); // C^P != \phi
};
//--------------------------------------------------------
QPatch::QPatch(){
	par=NULL;
	lev=0;
	icrs=0;	
	bndr=false; intr=false; extr=false;
	for(int i=0;i<4;i++) chld[i]=NULL;
};
int QPatch::isin(){
	if(intr) return(0);
	if(bndr) return(1);
	if(extr) return(2);
	return(-1);
};
QPatch::~QPatch(){
	//printf("QPatch destructed\n");
};
void QPatch::print(){
	px.print();
};
void QPatch::draw(char fn[128],char mode[3]){
	px.draw(fn,mode);
};
void QPatch::draw(FILE *fp){
	px.draw(fp);
};
void QPatch::draw(){
	px.draw();
};
void QPatch::set_lim(double xa[2], double xb[2]){
	px.set_Xa(xa[0],xa[1]);
	px.set_Xb(xb[0],xb[1]);
	px.setup();
};	
QPatch *new_QPatch(double Xa[2], double Wd[2]){
	QPatch *qp=(QPatch *)malloc(sizeof(QPatch));
	double Xb[2];
	Xb[0]=Xa[0]+Wd[0];
	Xb[1]=Xa[1]+Wd[1];
	qp->icrs=0;
	qp->set_lim(Xa,Xb);
	qp->bndr=false; qp->intr=false; qp->extr=false;
	for(int i=0;i<4;i++) qp->chld[i]=NULL;
	return(qp);
};
void new_QPatch(double Xa[2], double Wd[2],QPatch **qp_new){
	//QPatch *qp=(QPatch *)malloc(sizeof(QPatch));
	QPatch *qp;
	(*qp_new)=(QPatch *)malloc(sizeof(QPatch));
	qp=(*qp_new);
	double Xb[2];
	Xb[0]=Xa[0]+Wd[0];
	Xb[1]=Xa[1]+Wd[1];
	qp->icrs=0;
	qp->set_lim(Xa,Xb);
	qp->bndr=false; qp->intr=false; qp->extr=false;
	for(int i=0;i<4;i++) qp->chld[i]=NULL;
};
//--------------------------------------------------------
void translate_crs(int icrs){

	if(icrs==0) puts("A ^ B = phi");
	if(icrs==1) puts("A < B ");
	if(icrs==2) puts("B < A ");
	if(icrs==3) puts("A ^ B != phi");
};
void gather_leaves(QPatch *qp_par, int *count, QPatch *qp_leaves){
	if(qp_par->chld[0]==NULL){
		//qp_par->draw();
		//printf("idx=%d\n",*count);
		qp_leaves[*count]=(*qp_par);
		(*count)++;
	}else{
		for(int k=0;k<4;k++) gather_leaves(qp_par->chld[k],count,qp_leaves);
	};
};
int QtreeFind(QPatch *qp, double xf[2]){

	if(qp->chld[0]==NULL) return(qp->isin());

	double *Xa=qp->px.Xa;
	double *Wd=qp->px.Wd;
	int icell,indx[2];

	while(xf[0]-Xa[0]<0.0) xf[0]+=Wd[0];
	while(xf[1]-Xa[1]<0.0) xf[1]+=Wd[1];
	while(xf[0]-Xa[0]>Wd[0]) xf[0]-=Wd[0];
	while(xf[1]-Xa[1]>Wd[0]) xf[1]-=Wd[1];

	indx[0]=floor(2.*(xf[0]-Xa[0])/Wd[0]);
	indx[1]=floor(2.*(xf[1]-Xa[1])/Wd[1]);
	icell=indx[1]*2+indx[0];

	int iin=QtreeFind(qp->chld[icell],xf);
	return(iin);
};
int Qtree(QPatch *qp, Circ cr,int *count){

	int i,j,k;
	double Xa[2],Ya[2],Wd[2];

	int lev=qp->lev;

	int icrs=is_cross(qp->px, cr);
	if(lev > 6){
		//qp->draw();
		(*count)++;
		qp->icrs=icrs;
		return(lev);
	}

	//int icrs=is_cross(qp->px, cr);
	
	//printf("Relation of G to A is %d\n",icrs);
	//translate_crs(icrs);
	Wd[0]=qp->px.Wd[0]*0.5;
	Wd[1]=qp->px.Wd[1]*0.5;
	Xa[0]=qp->px.Xa[0];
	Xa[1]=qp->px.Xa[1];
	if(icrs>1){
		k=0;
		for(j=0; j<2; j++){
			Ya[1]=Xa[1]+Wd[1]*j;
		for(i=0; i<2; i++){
			Ya[0]=Xa[0]+Wd[0]*i;
			qp->chld[k]=new_QPatch(Ya,Wd);
			//new_QPatch(Ya,Wd,qp->chld[k]);
			qp->chld[k]->par=qp;
			qp->chld[k]->lev=lev+1;
			Qtree(qp->chld[k], cr, count);
			k++;
		}
		}
	}else{
		(*count)++;
		qp->icrs=icrs;
		//qp->draw();
	};
	return(qp->lev);
};
//---------------------------------------------------------------
double area(Ellip el1, Ellip el2, int lev_max, bool isect){

	Bbox bx;
	if(isect){
		bx=bbox_cross(el1.bbox,el2.bbox);
	}else{
		bx=bbox_union(el1.bbox,el2.bbox);
	};

	QPatch qp0;
	qp0.set_lim(bx.Xa,bx.Xb);
	int count=0;
	Qtree(&qp0,el1,el2,isect,&count,lev_max);
	QPatch *qp_leaves=(QPatch *)malloc(sizeof(QPatch)*count);
	count=0;
	gather_leaves(&qp0,&count,qp_leaves);

	int i,lev;
	double S=0.0;
	double ds0=bx.Wd[0]*bx.Wd[1];
	double *ds=(double *)malloc(sizeof(double)*(lev_max+1));
	ds[0]=ds0;
	for(i=1;i<=lev_max;i++){
		ds[i]=ds[i-1]*0.25;
	};
	for(i=0;i<count;i++){
		if(qp_leaves[i].extr) continue;
		lev=qp_leaves[i].lev;	
		if(qp_leaves[i].intr) S+=ds[lev];
		if(qp_leaves[i].bndr) S+=(0.5*ds[lev]);
	};
	free(qp_leaves);
	free(ds);
	clear_Qtree(&qp0);
	return(S);
};
void clear_Qtree(QPatch *qp){
	QPatch *par;
	if(qp->chld[0]!=NULL){
		for(int i=0;i<4;i++){
			clear_Qtree(qp->chld[i]);
			free(qp->chld[i]);
			qp->chld[i]=NULL;
		};
	};
	qp->intr=false;
	qp->bndr=false;
	qp->extr=false;
	qp->lev=0;
	qp->icrs=0;	
};
void clear_Qtree2(QPatch *qp){
	int i;
	if(qp->chld[0]!=NULL){
		for(i=0;i<4;i++) clear_Qtree(qp->chld[i]);
	};
	if(qp->lev!=0){
		free(qp);
	}else{
		qp->intr=false;
		qp->bndr=false;
		qp->extr=false;
		qp->icrs=0;	
		for(i=0;i<4;i++) qp->chld[i]=NULL;
	}
};

int Qtree(QPatch *qp, Ellip el1, Ellip el2, bool isect, int *count, int lev_max){

	int lev=qp->lev;
	Poly pl;
	pl.np=4;
	pl.xs=qp->px.xs;
	pl.ys=qp->px.ys;

	int ic;
	if(isect){ // set intersection
		qp->intr=true;
		qp->extr=false;
		ic=poly_cross(pl, el1);
		if(ic!=1) qp->intr=false; // intersection (AND)
		if(ic==0) qp->extr=true; 

		ic=poly_cross(pl, el2);
		if(ic!=1) qp->intr=false; // intersection (AND)
		if(ic==0) qp->extr=true; 
	}else{	// set sum  (union)
		qp->intr=false;
		qp->extr=true;
		ic=poly_cross(pl, el1);
		if(ic==1) qp->intr=true; // union (OR)
		if(ic!=0) qp->extr=false; 
		ic=poly_cross(pl, el2);
		if(ic==1) qp->intr=true; // union (OR)
		if(ic!=0) qp->extr=false; 
	};

	qp->bndr=false;
	if( !(qp->intr)){
	       if(!(qp->extr)){
		       qp->bndr=true;
	       }
	 }

	int icrs=2;
	if(qp->extr) icrs=0;
	if(qp->intr) icrs=1;
	if(qp->bndr) icrs=3;

	if(lev >= lev_max){
		(*count)++;
		qp->icrs=icrs;
		return(lev);
	}

	int i,j,k;
	double Xa[2],Ya[2],Wd[2];
	Wd[0]=qp->px.Wd[0]*0.5;
	Wd[1]=qp->px.Wd[1]*0.5;
	Xa[0]=qp->px.Xa[0];
	Xa[1]=qp->px.Xa[1];
	if(icrs>1){
		k=0;
		for(j=0; j<2; j++){
			Ya[1]=Xa[1]+Wd[1]*j;
		for(i=0; i<2; i++){
			Ya[0]=Xa[0]+Wd[0]*i;
			qp->chld[k]=new_QPatch(Ya,Wd);
			//new_QPatch(Ya,Wd,&(qp->chld[k]));
			qp->chld[k]->par=qp;
			qp->chld[k]->lev=lev+1;
			Qtree(qp->chld[k], el1,el2,isect, count,lev_max);
			k++;
		}
		}
	}else{
		(*count)++;
		qp->icrs=icrs;
		//qp->draw();
	};
	return(qp->lev);
};
int Qtree( 
	QPatch *qp, 	// quad-tree patch
	Ellip *els, int nelp, 	// ellipses 
	bool isect, 	// set operation true(=intersection)/false(=union)
	int *count, 	// number of leaves
	int lev_max,	// maximum tree height 
	Bbox unit_cell	// unit cell (to apply Periodic B.C.)
){

	int lev=qp->lev;
	Poly pl;
	pl.np=4;
	pl.xs=qp->px.xs;
	pl.ys=qp->px.ys;
	Ellip elp,el;

	int i1,i2,j1,j2;
	int p,i,j;
	double ux,uy;

	int ic;
	if(isect){ // set intersection
		qp->intr=true;
		qp->extr=false;
		for(p=0; p <nelp; p++){
			i1=0; i2=0;
			j1=0; j2=0;
			indx4PrdBC(&i1,&i2,&j1,&j2,els[p].bbox,unit_cell);
			for(i=i1;i<=i2;i++){
				ux=unit_cell.Wd[0]*i;
			for(j=j1;j<=j2;j++){
				uy=unit_cell.Wd[1]*j;
				elp=els[p];
				elp.slide(ux,uy);
				if(bbox_cross(elp.bbox,qp->px)){
					ic=poly_cross(pl, elp);
				}else{
					ic=0;
				};
				if(ic!=1) qp->intr=false; // intersection (AND)
				if(ic==0) qp->extr=true; 
			}
			}
		}

	}else{	// set sum  (union)
		qp->intr=false;
		qp->extr=true;
		for(p=0; p<nelp; p++){
			i1=0; i2=0;
			j1=0; j2=0;
			indx4PrdBC(&i1,&i2,&j1,&j2,els[p].bbox,unit_cell);
			for(i=i1;i<=i2;i++){
				ux=unit_cell.Wd[0]*i;
			for(j=j1;j<=j2;j++){
				uy=unit_cell.Wd[1]*j;
				elp=els[p];
				elp.slide(ux,uy);
				if(bbox_cross(elp.bbox,qp->px)){
					ic=poly_cross(pl, elp);
					if(ic==1) qp->intr=true; // union (OR)
					if(ic!=0) qp->extr=false; 
				}
			}
			}
		}
	};

	qp->bndr=false;
	if( !(qp->intr)){
	       if(!(qp->extr)){
		       qp->bndr=true;
	       }
	 }

	int icrs=2;
	if(qp->extr) icrs=0;
	if(qp->intr) icrs=1;
	if(qp->bndr) icrs=3;

	if(lev >= lev_max){
		(*count)++;
		qp->icrs=icrs;
		return(lev);
	}

	int k;
	double Xa[2],Ya[2],Wd[2];
	Wd[0]=qp->px.Wd[0]*0.5;
	Wd[1]=qp->px.Wd[1]*0.5;
	Xa[0]=qp->px.Xa[0];
	Xa[1]=qp->px.Xa[1];
	if(icrs>1){
		k=0;
		for(j=0; j<2; j++){
			Ya[1]=Xa[1]+Wd[1]*j;
		for(i=0; i<2; i++){
			Ya[0]=Xa[0]+Wd[0]*i;
			qp->chld[k]=new_QPatch(Ya,Wd);
			//new_QPatch(Ya,Wd,&(qp->chld[k]));
			qp->chld[k]->par=qp;
			qp->chld[k]->lev=lev+1;
			Qtree(qp->chld[k], els,nelp,isect, count,lev_max,unit_cell);
			k++;
		}
		}
	}else{
		(*count)++;
		qp->icrs=icrs;
	};
	return(qp->lev);
};
int Qtree(QPatch *qp, Ellip *els,int nelp, bool isect, int *count, int lev_max){

	int lev=qp->lev;
	Poly pl;
	pl.np=4;
	pl.xs=qp->px.xs;
	pl.ys=qp->px.ys;

//	printf("lev,lev_max=%d/%d\n",lev,lev_max);
	int ic,i;
	if(isect){ // set intersection
		qp->intr=true;
		qp->extr=false;
		for(i=0; i<nelp; i++){
			if(bbox_cross(els[i].bbox,qp->px)){
				ic=poly_cross(pl, els[i]);
			}else{
				ic=0;
			};
			if(ic!=1) qp->intr=false; // intersection (AND)
			if(ic==0) qp->extr=true; 
		}

	}else{	// set sum  (union)
		qp->intr=false;
		qp->extr=true;
		for(i=0; i<nelp; i++){
			if(bbox_cross(els[i].bbox,qp->px)){
				ic=poly_cross(pl, els[i]);
				if(ic==1) qp->intr=true; // union (OR)
				if(ic!=0) qp->extr=false; 
			}
		}
	};

	qp->bndr=false;
	if( !(qp->intr)){
	       if(!(qp->extr)){
		       qp->bndr=true;
	       }
	 }

	int icrs=2;
	if(qp->extr) icrs=0;
	if(qp->intr) icrs=1;
	if(qp->bndr) icrs=3;

	if(lev >= lev_max){
		(*count)++;
		qp->icrs=icrs;
		return(lev);
	}

	int j,k;
	double Xa[2],Ya[2],Wd[2];
	Wd[0]=qp->px.Wd[0]*0.5;
	Wd[1]=qp->px.Wd[1]*0.5;
	Xa[0]=qp->px.Xa[0];
	Xa[1]=qp->px.Xa[1];
	if(icrs>1){
		k=0;
		for(j=0; j<2; j++){
			Ya[1]=Xa[1]+Wd[1]*j;
		for(i=0; i<2; i++){
			Ya[0]=Xa[0]+Wd[0]*i;
			qp->chld[k]=new_QPatch(Ya,Wd);
			//new_QPatch(Ya,Wd,&(qp->chld[k]));
			qp->chld[k]->par=qp;
			qp->chld[k]->lev=lev+1;
			Qtree(qp->chld[k], els,nelp,isect, count,lev_max);
			k++;
		}
		}
	}else{
		(*count)++;
		qp->icrs=icrs;
	};
	return(qp->lev);
};
//---------------------------------------------------------------
int Qtree(QPatch *qp, Solid sld,int *count, int lev_max){

	int i,j,k,icrs;
	double Xa[2],Ya[2],Wd[2];

	int lev=qp->lev;
	Poly pl;
	pl.np=4;
	pl.xs=qp->px.xs;
	pl.ys=qp->px.ys;
	int ic;
	bool init=false;

	// Tile union of all sets (particles);
	int isum=0;
	for(i=0;i<sld.nelp;i++){
		if(sld.isect[i]) continue;
		if(isum==0){
			qp->intr=false;
			qp->extr=true;
			init=true;
		}
		if(bbox_cross(sld.els[i].bbox,qp->px)){
			ic=poly_cross(pl, sld.els[i]);
			if(ic==1) qp->intr=true; // union (OR)
			if(ic!=0) qp->extr=false; 
		}
		isum++;
	}
	// Tile intersection of all sets (particles)
	if(!init){
		qp->intr=true;
		qp->extr=false;
	}
	for(i=0;i<sld.nelp;i++){
	if(sld.isect[i]){
		if(bbox_cross(sld.els[i].bbox,qp->px)){
			ic=poly_cross(pl, sld.els[i]);
		}else{
			ic=0;
		}
			if(ic!=1) qp->intr=false; // intersection (AND)
			if(ic==0) qp->extr=true; 
	}
	}

	qp->bndr=false;
	if( !(qp->intr)){
	       if(!(qp->extr)){
		       qp->bndr=true;
	       }
	 }

	icrs=2;
	if(qp->extr) icrs=0;
	if(qp->intr) icrs=1;
	if(qp->bndr) icrs=3;

	if(lev >= lev_max){
		(*count)++;
		qp->icrs=icrs;
		return(lev);
	}

	//translate_crs(icrs);
	Wd[0]=qp->px.Wd[0]*0.5;
	Wd[1]=qp->px.Wd[1]*0.5;
	Xa[0]=qp->px.Xa[0];
	Xa[1]=qp->px.Xa[1];
	if(icrs>1){
		k=0;
		for(j=0; j<2; j++){
			Ya[1]=Xa[1]+Wd[1]*j;
		for(i=0; i<2; i++){
			Ya[0]=Xa[0]+Wd[0]*i;
			//new_QPatch(Ya,Wd,&(qp->chld[k]));
			qp->chld[k]=new_QPatch(Ya,Wd);
			qp->chld[k]->par=qp;
			qp->chld[k]->lev=lev+1;
			Qtree(qp->chld[k], sld, count,lev_max);
			k++;
		}
		}
	}else{
		(*count)++;
		qp->icrs=icrs;
		//qp->draw();
	};
	return(qp->lev);
};
//---------------------------------------------------------
void Tree4::setup(Ellip el1, Ellip el2, bool set_opr, int LevMax){

	int count;
	lev_max=LevMax;
	isect=set_opr;

	count=0;
	el1.set_bbox();
	el2.set_bbox();
	Bbox bx;
	if(isect){
		bx=bbox_cross(el1.bbox,el2.bbox);
	}else{
		bx=bbox_union(el1.bbox,el2.bbox);
	};

	qp0.set_lim(bx.Xa, bx.Xb);
	Qtree(&qp0,el1,el2,isect,&count,lev_max);
	n_leaves=count;

	leaves=(QPatch *)malloc(sizeof(QPatch)*n_leaves);
	count=0;
	gather_leaves(&qp0,&count,leaves);
	ready=true;
};
void Tree4::setup(Ellip *els, int nelp, bool set_opr, int LevMax){

	int i,count;
	lev_max=LevMax;
	count=0;
	isect=set_opr;


	for(i=0;i<nelp;i++) els[i].set_bbox();
	Bbox bx=els[0].bbox;
	if(isect){
		for(i=1;i<nelp;i++) bx=bbox_cross(bx,els[i].bbox);
	}else{
		for(i=1;i<nelp;i++) bx=bbox_union(bx,els[i].bbox);
	}

	qp0.set_lim(bx.Xa, bx.Xb);
	Qtree(&qp0,els,nelp,isect,&count,lev_max);
	n_leaves=count;

	leaves=(QPatch *)malloc(sizeof(QPatch)*n_leaves);
	count=0;
	gather_leaves(&qp0,&count,leaves);
	ready=true;

}
void Tree4::setup(Ellip *els, int nelp, bool set_opr, int LevMax, Bbox bx){

	int i,count;
	lev_max=LevMax;
	count=0;
	isect=set_opr;

	qp0.set_lim(bx.Xa, bx.Xb);
	Qtree(&qp0,els,nelp,isect,&count,lev_max,bx);
	n_leaves=count;

	leaves=(QPatch *)malloc(sizeof(QPatch)*n_leaves);
	count=0;
	gather_leaves(&qp0,&count,leaves);
	ready=true;

}
void Tree4::setup(Solid sld,int LevMax){

	int count;
	lev_max=LevMax;
	count=0;

	qp0.set_lim(sld.bbox.Xa,sld.bbox.Xb);
	Qtree(&qp0,sld,&count,lev_max);
	n_leaves=count;

	leaves=(QPatch *)malloc(sizeof(QPatch)*n_leaves);
	count=0;
	gather_leaves(&qp0,&count,leaves);
	ready=true;

};
Tree4::Tree4(){
	lev_max=0;
	isect=true;
	ready=false;
};
void Tree4::draw(){
	int i;
	char mode[3]="a";
	char fni[128],fnb[128],fne[128];
	FILE *fp1,*fp2,*fp3;

	if(!ready){
		printf("Tree strcuture has yet to been build");
	}else{
		sprintf(fni,"qtree_in.out");
		sprintf(fnb,"qtree_bnd.out");
		sprintf(fne,"qtree_ex.out");
		fp1=fopen(fni,"w");
		fp2=fopen(fnb,"w");
		fp3=fopen(fne,"w");
		for(i=0;i<n_leaves;i++){
			if(leaves[i].intr) leaves[i].draw(fp1);
			if(leaves[i].bndr) leaves[i].draw(fp2);
			if(leaves[i].extr) leaves[i].draw(fp3);
		}
		fclose(fp1);
		fclose(fp2);
		fclose(fp3);
	};
};
void Tree4::count(){
	int i;
	nint=0;
	next=0;
	nbnd=0;
	if(!ready){
		printf("Tree strcuture has yet to been build");
	}else{
		for(i=0;i<n_leaves;i++){
			if(leaves[i].intr) nint++;
			if(leaves[i].bndr) nbnd++;
			if(leaves[i].extr) next++;
		}
	};
	printf("nint=%d, nbnd=%d, next=%d\n",nint,nbnd,next);
};
void Tree4::set_grid_params(){
	Nx=pow(2,lev_max);
	Ny=Nx;
	for(int i=0;i<2;i++){
		Xa[i]=qp0.px.Xa[i];
		Xb[i]=qp0.px.Xb[i];
		Wd[i]=qp0.px.Wd[i];
	}
	dx[0]=Wd[0]/Nx;
	dx[1]=Wd[1]/Ny;
	//Nx++; Ny++;
};
int Tree4::grid_type(int id, int jd, int cnct[4]){
	// 0: interior grid
	// 1: boundary grid
	// 2: exterior grid
	// -1: error
	
	Tree4::set_grid_params();
	while(id<0) id+=Nx; 
	while(jd<0) jd+=Ny; 
	while(id>=Nx) id-=Nx; 
	while(jd>=Ny) jd-=Ny; 

	double xf[2],x0[2];
	x0[0]=Xa[0]+dx[0]*id;
	x0[1]=Xa[1]+dx[1]*jd;

	int i,j,ityp;
	int ntyp[3]={0,0,0};
	int isgn[4]={-1, 1, 1,-1};
	int jsgn[4]={-1,-1, 1, 1};
	int isum=0;

	for(i=0;i<4;i++) cnct[i]=0;
	for(j=0;j<2;j++){
		xf[1]=x0[1]+jsgn[isum]*dx[1]*0.5;
	for(i=0;i<2;i++){
		xf[0]=x0[0]+isgn[isum]*dx[0]*0.5;
		ityp=QtreeFind(&qp0,xf);
		if(ityp!=-1) ntyp[ityp]++;
		if(ityp==1){
			cnct[isum]=1;
			cnct[(isum+1)%4]=1;
		};
		isum++;
	}
	};

	if(ntyp[1]>0){
//		printf("%lf %lf\n",x0[0],x0[1]);
		return(1);
	};
	if(ntyp[0]>ntyp[2]) return(0);
	if(ntyp[2]>ntyp[1]) return(2);
	return(-1);
};
double Tree4::area(){
	int i,lev;
	double S=0.0;
	double ds0=qp0.px.area();
	double *ds=(double *)malloc(sizeof(double)*(lev_max+1));
	ds[0]=ds0;
	for(i=1;i<=lev_max;i++) ds[i]=ds[i-1]*0.25;
	for(i=0;i<n_leaves;i++){
		if(leaves[i].extr) continue;
		lev=leaves[i].lev;	
		if(leaves[i].intr) S+=ds[lev];
		if(leaves[i].bndr) S+=(0.5*ds[lev]);
	};
	free(ds);
	return(S);
};
void Tree4::clean(){
	if(ready){
		free(leaves);
		clear_Qtree(&qp0);
		//clear_Qtree2(&qp0);
		ready=false;
	};
};
//---------------------------------------------------------
Ellip::Ellip(){
	xc[0]=0.0; xc[1]=0.0;
	radi[0]=1.0; radi[1]=1.0;
	phi=0.0;
};
double Ellip::area(){
	return(radi[0]*radi[1]*4.0*atan(1.0));
};
void Ellip::set_bbox(){
	double cosp=cos(phi);
	double sinp=sin(phi);
	double a2=radi[0]*radi[0];
	double b2=radi[1]*radi[1];

	cosp*=cosp;
	sinp*=sinp;
	double xmax=sqrt(a2*cosp+b2*sinp);
	double ymax=sqrt(a2*sinp+b2*cosp);

	bbox.set_Xa(xc[0]-xmax,xc[1]-ymax);
	bbox.set_Xb(xc[0]+xmax,xc[1]+ymax);

};
void Ellip::scale(double s){
	radi[0]*=s;
	radi[1]*=s;
};
void Ellip::slide(double ux, double uy, Bbox unit_cell){
	xc[0]+=ux;
	xc[1]+=uy;

	double *Xa=unit_cell.Xa;
	double *Xb=unit_cell.Xb;
	double *Wd=unit_cell.Wd;

	while( xc[0] < Xa[0]) xc[0]+=Wd[0];
	while( xc[1] < Xa[1]) xc[1]+=Wd[1];
	while( xc[0] > Xb[0]) xc[0]-=Wd[0];
	while( xc[1] > Xb[1]) xc[1]-=Wd[1];

	Ellip::set_bbox();
};
void Ellip::slide(double ux, double uy){
	xc[0]+=ux;
	xc[1]+=uy;
	bbox.slide(ux,uy);
	//Ellip::set_bbox();
};
void Ellip::draw(int np){
	double dth=8.0*atan(1.0)/np;
	double th;
	double x,y,X,Y;
	double cosp=cos(phi),sinp=sin(phi);
	for(int i=0;i<=np;i++){
		th=i*dth;
		x=radi[0]*cos(th);
		y=radi[1]*sin(th);
		X= cosp*x-sinp*y+xc[0];
		Y= sinp*x+cosp*y+xc[1];
	       	printf("%lf %lf\n",X,Y);
	}
	printf("\n");
};
void Ellip::draw(char fn[128], int np, char mode[3]){
	FILE *fp=fopen(fn,mode);
	double dth=8.0*atan(1.0)/np;
	double th;
	double x,y,X,Y;
	double cosp=cos(phi),sinp=sin(phi);
	for(int i=0;i<=np;i++){
		th=i*dth;
		x=radi[0]*cos(th);
		y=radi[1]*sin(th);
		X= cosp*x-sinp*y+xc[0];
		Y= sinp*x+cosp*y+xc[1];
	       	fprintf(fp,"%lf %lf\n",X,Y);
	}
	fprintf(fp,"\n");

	fclose(fp);
};
void Ellip::set_xc(double x, double y){
	xc[0]=x;
	xc[1]=y;
};
void Ellip::set_radi(double r1, double r2){
	radi[0]=r1;
	radi[1]=r2;
};
void Ellip::set_phi(double ang){
	double pi=4.0*atan(1.0);
	phi=ang/180.*pi;
};
bool Ellip::is_in(double xf[2]){
	double yf[2],Y[2];

	yf[0]=xf[0]-xc[0];
	yf[1]=xf[1]-xc[1];

	double cosp=cos(phi),sinp=sin(phi);
	Y[0]=cosp*yf[0]+sinp*yf[1];
	Y[1]=-sinp*yf[0]+cosp*yf[1];
	Y[0]/=radi[0];
	Y[1]/=radi[1];

	bool iin=true;

	if(Y[0]*Y[0]+Y[1]*Y[1]>1.0) iin=false;
	return(iin);
};
Poly::Poly(){
	alocd=false;
};
Poly::Poly(int n){
	np=n;
	Poly::mem_alloc();
	xg[0]=0.0;
	xg[1]=0.0;
};
Poly::~Poly(){
//	Thie mem_free destructor does not work
//		for a reason not known !! 
//	if(alocd){
//	      free(xs);
//	      free(ys);
//	}
};
void Poly::mem_alloc(){
	xs=(double *)malloc(sizeof(double)*np);
	ys=(double *)malloc(sizeof(double)*np);
	xg[0]=0.0;
	xg[1]=0.0;
	alocd=true;
};
void Poly::mem_free(){
	free(xs);
	free(ys);
	alocd=false;
};
void Poly::set_center(){
	for(int i=0;i<np;i++){
		xg[0]+=xs[i];
		xg[1]+=ys[i];
	}
	xg[0]/=np;
	xg[1]/=np;
}
void Poly::set_bbox(){
	double xmax=xs[0];
	double xmin=xs[0];
	double ymax=ys[0];
	double ymin=ys[0];
	for(int i=1;i<np;i++){
		if(xs[i]>xmax) xmax=xs[i];
		if(xs[i]<xmin) xmin=xs[i];

		if(ys[i]>ymax) ymax=ys[i];
		if(ys[i]<ymin) ymin=ys[i];
	};

	bbox.set_Xa(xmin,ymin);
	bbox.set_Xb(xmax,ymax);
};
void Poly::draw(){
	for(int i=0;i<=np;i++){
		printf("%lf %lf\n",xs[i%np],ys[i%np]);
	};
	printf("\n");
};
void Poly::draw(char fn[128],char mode[3]){
	FILE *fp=fopen(fn,mode);
	for(int i=0;i<=np;i++){
		fprintf(fp,"%lf %lf\n",xs[i%np],ys[i%np]);
	};
	fprintf(fp,"\n");
	fclose(fp);
};
bool Poly::is_in(double xf[2]){
	int i;
	Vec2 x1,x2,yf,tx,nx;
	yf.set(xf[0],xf[1]);
	for(i=0;i<np;i++){
		x1.set(xs[i],ys[i]);	
		x2.set(xs[(i+1)%np],ys[(i+1)%np]);	
		tx=vdiff(x2,x1);
		nx.set(tx.x[1],-tx.x[0]);
		if(iprod(nx, vdiff(yf,x1))>0.0){
			return(false);
		}
	};
	return(true);
};
void Poly::slide(double ux, double uy){
	for(int i=0;i<np;i++){
		xs[i]+=ux;
		ys[i]+=uy;
	}
};
void Poly::rotate(int node_num, double th){

	Vec2 xc,ri;
	if(node_num<0){
		xc.set(xg);	// centroide
	}else{
		xc.set(xs[node_num%np],ys[node_num%np]);
	};
		

	for(int i=0; i<np; i++){
		ri.set(xs[i],ys[i]);
		ri=vdiff(ri,xc);
		ri.rotate(th);
		ri=vsum(ri,xc);
		xs[i]=ri.x[0];
		ys[i]=ri.x[1];
	}
};
void Poly::rotate(double xr0[2], double th){
	Vec2 xc(xr0);
	Vec2 ri;

	for(int i=0; i<np; i++){
		ri.set(xs[i],ys[i]);
		ri=vdiff(ri,xc);
		ri.rotate(th);
		ri=vsum(ri,xc);
		xs[i]=ri.x[0];
		ys[i]=ri.x[1];
	}
};

int poly_cross(Poly A, Poly B){
	int i,j;
	double xf[2];

	int intr=0;
	int extr=0;
	for(i=0;i<A.np;i++){
		xf[0]=A.xs[i];
		xf[1]=A.ys[i];
		if(B.is_in(xf)) intr++;
	};

	if(intr==A.np) return(1); // A < B

	int jntr=0;
	for(i=0;i<B.np;i++){
		xf[0]=B.xs[i];
		xf[1]=B.ys[i];
		if(A.is_in(xf)) jntr++;
	};

	if(jntr==B.np) return(2); // B < A
	if((intr+jntr)==0) return(0); // A^B=\phi

	return(3); // A^B != \phi
};
int poly_cross(Poly A, Circ B){
	int intr=0,extr=0;
	int i;
	double x1[2],x2[2];
	double rmin,rmax;
	int np=A.np;

	for(i=0; i<np; i++){
		x1[0]=A.xs[i]; x1[1]=A.ys[i];
		x2[0]=A.xs[(i+1)%np]; x2[1]=A.ys[(i+1)%np];
		distP2L(B.xc,x1,x2,&rmin,&rmax);	
		if(rmax < B.radi) intr++;
		if(rmin > B.radi) extr++;
	}
	if(intr==np) return(1); // A< B
	if(extr==np){
		if(A.is_in(B.xc)){
			return(2); //B <A
		}else{
			return(0); // A ^ B = \phi
		}
	}	
	return(3); // A ^ B != \phi
};
int poly_cross(Poly A, Ellip B){

	int np=A.np;
	double phi=B.phi;
	double cosp=cos(phi),sinp=sin(phi);

	Poly Ad(np);
	Circ Bd;
	Bd.xc[0]=0.0;
	Bd.xc[1]=0.0;
	Bd.radi=1.0;
	//Bd.draw("temp.dat",100,"a");

	double *xcod=Ad.xs;
	double *ycod=Ad.ys;

	double x,y;
	for(int i=0;i<np;i++){
		x=A.xs[i]-B.xc[0];
		y=A.ys[i]-B.xc[1];
		xcod[i]= cosp*x + sinp*y;
		ycod[i]=-sinp*x + cosp*y;
		xcod[i]/=B.radi[0];
		ycod[i]/=B.radi[1];
	};

	//Ad.draw("temp.dat","a");
	int icrs=poly_cross(Ad,Bd);
	Ad.mem_free();
	return(icrs);
};
Solid::Solid(){
	Solid::init_rand(-1);
};
Solid::Solid(int n){
	nelp=n;
	els=(Ellip *)malloc(sizeof(Ellip)*nelp);
	isect=(bool *)malloc(sizeof(bool)*nelp);
	Solid::init_rand(-1);
};
Solid::Solid(int n, double Wd[2]){

	nelp=n;
	els=(Ellip *)malloc(sizeof(Ellip)*nelp);
	isect=(bool *)malloc(sizeof(bool)*nelp);
	Solid::init_rand(-1);

	double Xa[2],Xb[2];
	Xa[0]=0.0; Xa[1]=0.0;	
	Xb[0]=Xa[0]+Wd[0];
	Xb[1]=Xa[1]+Wd[1];
	bbox.setup(Xa,Xb);

	double x,y,phi;
	double PI=4.0*atan(1.0);
	double r0=Wd[0]*0.03,aspect=0.6;

	double S0=0.0;
	for(int i=0;i<nelp;i++){
		x=Urnd(mt)*Wd[0]+Xa[0];
		y=Urnd(mt)*Wd[1]+Xa[1];
		els[i].set_xc(x,y);
		els[i].phi=Urnd(mt)*PI;
		els[i].set_radi(r0,r0*aspect);
		isect[i]=false;
		els[i].set_bbox();
		S0+=els[i].area();
	};
	printf("Number of particles=%d\n",nelp);
	double psi=S0/bbox.area();
	printf("Max. packing density=%lf%%\n",psi*100.);
	printf("(min. porosity=%lf\n",1-psi);
};
void Solid::init_rand(int seed){
	mt=std::mt19937_64(seed);
	Urnd=std::uniform_real_distribution<double>(0.0,1.0); // Uniform distribution(min,max)
	Grnd=std::normal_distribution<double>(0.0,1.0); // Normal distribution(mean,stdev)	
};
void Solid::draw(char fn[128],int ndat){
	char md[3];
	for(int i=0;i<nelp;i++){
		sprintf(md,"a");
		if(i==0) sprintf(md,"w");
		els[i].draw(fn,ndat,md);
	}
};
double Solid::MC(Temp_Hist TH){
	double PI=4.0*atan(1.0);
	int ip;
	double ux,uy,dphi;
	double dE,dE_sum=0.0;
	double prb,alph,beta;
	double tau=TH.tau();
	double a1=1.0, a2=0.1;
	alph=(1.-tau)*a1+tau*a2;	// a1 --> a2
	//beta=(1.-tau)*1.0+tau*0.5;	// 1.0 --> 0.5 
	beta=1.0;
	for(ip=0; ip<nelp; ip++){
		if(ip%100==0) printf("ip=%d\n",ip);
		ux=0.0; uy=0.0; dphi=0.0;
		if(ip %2==0){
			ux=(2.*(Urnd(mt)-0.5)*bbox.Wd[0]*0.5 )*alph;
			uy=(2.*(Urnd(mt)-0.5)*bbox.Wd[1]*0.5 )*alph;
		}else{
			dphi=PI*Urnd(mt)*beta;
		};
	       	//dE=perturb(ip,ux,uy,dphi);
	       	dE=perturb_periodic(ip,ux,uy,dphi);
		prb=exp(-dE/TH.Temp);
		if(prb>1.0) prb=1.0;
		if(Urnd(mt) <=prb){	//accept
			els[ip].phi+=dphi;
			els[ip].slide(ux,uy,bbox); 
			dE_sum+=dE;
		};
	}
	return(dE_sum);
};
double Solid::perturb_periodic(int p, double ux, double uy, double dphi){

	Ellip elp,elq,elr,elr0;
	Tree4 tr4;
	int q,lev_max=6;
	int i,j,k,l;
	int i1,i2,j1,j2;
	int k1,k2,l1,l2;
	double Ux,Uy,Vx,Vy;
	bool pq;
	double dEp=0.0,dEm=0.0,S=0.0;

	i1=0; i2=0;
	j1=0; j2=0;
	indx4PrdBC(&i1,&i2,&j1,&j2,els[p].bbox,bbox);

	for(i=i1;i<=i2;i++){
		Ux=bbox.Wd[0]*i;
	for(j=j1;j<=j2;j++){
		Uy=bbox.Wd[1]*j;
		elp=els[p];
		elp.slide(Ux,Uy);
		for(q=0; q<nelp; q++){
			k1=0; k2=0;
			l1=0; l2=0;
			indx4PrdBC(&k1,&k2,&l1,&l2,els[q].bbox,bbox);
			for(k=k1;k<=k2;k++){
				Vx=bbox.Wd[0]*k;
			for(l=l1;l<=l2;l++){
				if(abs(p-q)+abs(k-i)+abs(j-l)==0) continue;
				Vy=bbox.Wd[1]*l;
				elq=els[q];
				elq.slide(Vx,Vy); 
				pq=bbox_cross(elp,elq);
				if(pq){
					tr4.setup(elp,elq,true,lev_max);
					S=tr4.area();
					tr4.clean();
					dEm+=S;
				}
			}	// end_l
			}	// end_k
		} // end_q
	}	// end_j
	}	// end_i


	S=0.0;
	elr=els[p]; // copy p-th ellipse before perturbation
	elr.phi+=dphi;
	elr.slide(ux,uy,bbox); 

	i1=0; i2=0;
	j1=0; j2=0;
	indx4PrdBC(&i1,&i2,&j1,&j2,elr.bbox,bbox);
	for(i=i1;i<=i2;i++){
		Ux=bbox.Wd[0]*i;
	for(j=j1;j<=j2;j++){
		Uy=bbox.Wd[1]*j;
		elp=elr;
		elp.slide(Ux,Uy);
		for(q=0; q<nelp; q++){
			k1=0; k2=0;
			l1=0; l2=0;
			indx4PrdBC(&k1,&k2,&l1,&l2,els[q].bbox,bbox);
			for(k=k1;k<=k2;k++){
				Vx=bbox.Wd[0]*k;
			for(l=l1;l<=l2;l++){
				if(abs(p-q)+abs(k-i)+abs(j-l)==0) continue;
				Vy=bbox.Wd[1]*l;
				elq=els[q];
				elq.slide(Vx,Vy); 
				pq=bbox_cross(elp,elq);
				if(pq){
					tr4.setup(elp,elq,true,lev_max);
					S=tr4.area();
					tr4.clean();
					dEp+=S;
				}
			}	// end_l
			}	// end_k
		} // end_q
	}	// end_j
	}	// end_i

	//printf("dE=%lf (dEp,dEm)=%lf %lf\n",dEp-dEm,dEp,dEm);
	return(dEp-dEm);

}
double Solid::perturb(int p, double ux, double uy, double dphi){

	int q,lev_max=6;
	bool pq;
	Ellip elp=els[p];
	Ellip elq,elr;
	double dEp=0.0,dEm=0.0,S=0.0;
	Tree4 tr4;

	for(q=0; q<nelp; q++){
		if(q==p) continue;
		pq=bbox_cross(elp,els[q]);
		if(pq){
			elq=els[q];
			tr4.setup(elp,elq,true,lev_max);
			S=tr4.area();
			tr4.clean();
			dEm+=S;
		}
	}	
	elr=elp;	// copy p-th ellipse before perturbation
	elr.phi+=dphi;
	elr.slide(ux,uy,bbox); 
	//fflush(stdout);

	S=0.0;
	for(q=0; q<nelp; q++){
		if(q==p) continue;
		pq=bbox_cross(elr,els[q]);
		if(pq){
			elq=els[q];
			tr4.setup(elr,elq,true,lev_max);
			S=tr4.area();
			tr4.clean();
			dEp+=S;
		}
	}	
	//printf("dE=%lf (dEp,dEm)=%lf %lf\n",dEp-dEm,dEp,dEm);
	return(dEp-dEm);
};
double Solid::area(int lev_max){
	Tree4 tr4;
	//tr4.setup(els,nelp,false,lev_max);
	tr4.setup(els,nelp,false,lev_max,bbox);
	double S=tr4.area();
	tr4.clean();

	double si=0.0;
	int lev=lev_max-3;
	if(lev <=1) lev=2;
	for(int i=0;i<nelp;i++){
		tr4.setup(els+i,1,false,lev);
		si+=tr4.area();
		tr4.clean();
	};
	printf("Overlap =%lf\n",si-S);
	double phi=S/bbox.area();
	printf("Packng Density =%lf%%\n",S*100.);
	printf("Porosity sity =%lf%%\n",1-S);
	return(S);
};

void indx4PrdBC(		// retrun indices to apply periodic B.C.
	int *i1, int *i2, 	// 1st index (start,end)
	int *j1, int *j2, 	// 2nd index (start,end)
	Bbox bx, 		// bounding box of geoemetric object
	Bbox B0			// Unit cell
){
	*i1=0; *i2=0;
	*j1=0; *j2=0;

	if(bx.Xa[0] < B0.Xa[0]) *i2= 1;
	if(bx.Xb[0] > B0.Xb[0]) *i1=-1;
	if(bx.Xa[1] < B0.Xa[1]) *j2= 1;
	if(bx.Xb[1] > B0.Xb[1]) *j1=-1;
};
#if DB==6	//  introducing periodic BC
int main(){
	double Xa[2]={0.0,0.0};
	double Xb[2]={1.0,1.0};
	Bbox BX,bx;
	BX.setup(Xa,Xb);

	Ellip el,elq;

	el.set_xc(1.1,0.3);
	el.set_radi(0.2,0.1);
	el.set_phi(30.0);
	el.set_bbox();

	int i,j;
	int i1,i2;
	int j1,j2;
	double ux,uy;

	indx4PrdBC(&i1,&i2,&j1,&j2,el.bbox,BX);
	//printf("(i1,i2)=(%d,%d)\n",i1,i2);
	//printf("(j1,j2)=(%d,%d)\n",j1,j2);
	for(i=i1;i<=i2;i++){
		ux=BX.Wd[0]*i;	
	for(j=j1;j<=j2;j++){
		uy=BX.Wd[1]*j;	
		elq=el;
		elq.slide(ux,uy);
		elq.draw(50);
		elq.bbox.draw();
	}
	};

	el.set_xc(1.1,-0.1);
	el.set_phi(-30.0);
	el.set_bbox();
	indx4PrdBC(&i1,&i2,&j1,&j2,el.bbox,BX);
	//printf("(i1,i2)=(%d,%d)\n",i1,i2);
	//printf("(j1,j2)=(%d,%d)\n",j1,j2);
	for(i=i1;i<=i2;i++){
		ux=BX.Wd[0]*i;	
	for(j=j1;j<=j2;j++){
		uy=BX.Wd[1]*j;	
		elq=el;
		elq.slide(ux,uy);
		elq.draw(50);
		elq.bbox.draw();
	}
	}

	BX.draw();

	return(0);
};
#endif
#if DB==5	//  testing bounding box
int main(){
	Ellip el;
	el.set_xc(1.0,0.5);
	el.set_radi(2.0,0.5);
	el.scale(0.5);
	el.set_phi(130.0);
	el.set_bbox();
	el.bbox.draw();
	el.draw(50);

	Poly plA(5);
	plA.xs[0]=0.0; plA.ys[0]=0.0;
	plA.xs[1]=1.0; plA.ys[1]=0.0;
	plA.xs[2]=1.0; plA.ys[2]=1.0;
	plA.xs[3]=0.5; plA.ys[3]=1.5;
	plA.xs[4]=0.0; plA.ys[4]=0.8;
	plA.rotate(1,30.0);
	plA.slide(2.0,2.0);
	plA.draw();
	plA.set_bbox();
	plA.bbox.draw();
	Bbox bbu=bbox_union(el.bbox,plA.bbox);
	Bbox bbc=bbox_cross(el.bbox,plA.bbox);

	return(0);
};
#endif
#if DB==4
int main(){
	char fname[128]="temp.dat";
	char md[3]="w";

	Poly plA(5);
	plA.xs[0]=0.0; plA.ys[0]=0.0;
	plA.xs[1]=1.0; plA.ys[1]=0.0;
	plA.xs[2]=1.0; plA.ys[2]=1.0;
	plA.xs[3]=0.5; plA.ys[3]=1.5;
	plA.xs[4]=0.0; plA.ys[4]=0.8;
	plA.set_center();
	plA.draw(fname,md);

	Ellip elB;
	elB.set_xc(-0.5,-0.5);
	elB.set_radi(1.5,0.5);


	sprintf(md,"a");
	double phi=45.0;
	for(int i=0;i<10;i++){
		elB.set_phi(phi);
		printf("icrs=%d\n",poly_cross(plA,elB));
		elB.draw(fname,100,md);
		//phi+=15.0;
		elB.radi[0]*=1.1;
		elB.radi[1]*=1.1;
	};

	return(0);
};
#endif
#if DB==3
int main(){
	char fname[128]="temp.dat";
	char md[3]="w";

	Poly plA(5);
	plA.xs[0]=0.0; plA.ys[0]=0.0;
	plA.xs[1]=1.0; plA.ys[1]=0.0;
	plA.xs[2]=1.0; plA.ys[2]=1.0;
	plA.xs[3]=0.5; plA.ys[3]=1.5;
	plA.xs[4]=0.0; plA.ys[4]=0.8;
	plA.set_center();

	Circ crB;
	crB.xc[0]=0.0;
	crB.xc[1]=0.0;
	crB.radi=0.25;
	plA.draw(fname,md);
	sprintf(md,"a");
	for(int i=0;i<10;i++){
		crB.draw(fname,100,md);
		printf("icrs=%d\n",poly_cross(plA,crB));
		//crB.radi+=0.1;
		crB.xc[0]+=0.25;
		crB.xc[1]+=0.25;
	};
	return(0);
};
#endif
#if DB==2
int main(){
	Poly plA(5);
	Poly plB(4);

	plA.xs[0]=0.0; plA.ys[0]=0.0;
	plA.xs[1]=1.0; plA.ys[1]=0.0;
	plA.xs[2]=1.0; plA.ys[2]=1.0;
	plA.xs[3]=0.5; plA.ys[3]=1.5;
	plA.xs[4]=0.0; plA.ys[4]=0.8;
	plA.set_center();


	plB.xs[0]=0.5; plB.ys[0]=0.0;
	plB.xs[1]=0.0; plB.ys[1]=0.5;
	plB.xs[2]=-0.5; plB.ys[2]=0.0;
	plB.xs[3]=-0.5; plB.ys[3]=-0.5;
	plB.set_center();

	double xf[2];

	char fname[128]="temp.dat";
	char md[3]="w";
	plA.draw(fname,md);
	plB.slide(-0.3,-0.3);

	sprintf(md,"a");
	double pi=4.*atan(1.0);
	double dth=10.0/180.*pi;
	xf[0]=0.0;
	xf[1]=0.0;
	for(int i=0;i<10;i++){
		//plB.slide(0.1,0.1);

		plA.draw(fname,md);
		plB.draw(fname,md);

		printf("icrs=%d\n",poly_cross(plA,plB));
		//plA.rotate(xf,dth);
		plA.rotate(1,dth);
	}
	return(0);
};
#endif
#if DB==1
int main(){
	Ellip el;
	el.set_xc(1.0,0.5);
	el.set_radi(2.0,0.5);
	el.set_phi(130.0);
	char fn[128]="temp.dat";
	char md[3]="w";
	el.draw(fn,100,md);

	double xf[2];
	double dx=0.1;
	double pi=4.0*atan(1.0);
	double th=15.0/180.*pi;
	for(int i=0;i<50;i++){
		xf[0]=(dx*i)*cos(th);	
		xf[1]=(dx*i)*sin(th);	
		printf("%lf %lf %d\n",xf[0],xf[1],el.is_in(xf));
	};	
	return(0);
};
#endif
