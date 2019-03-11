#define DB 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vec2.h"
#include "set2d.h"

using namespace std;
void Bbox::set_Xa(double x, double y){
	Xa[0]=x; Xa[1]=y;
};
void Bbox::set_Xb(double x, double y){
	Xb[0]=x; Xb[1]=y;
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
	return(b3);
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
void QPatch::print(){
	px.print();
};
void QPatch::draw(char fn[128],char mode[3]){
	px.draw(fn,mode);
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
int Qtree(QPatch *qp, Solid sld,int *count){

	int i,j,k,icrs;
	double Xa[2],Ya[2],Wd[2];

	int lev=qp->lev;

	icrs=0;
	Ellip elp;
	Poly pl;
	pl.np=4;
	pl.xs=qp->px.xs;
	pl.ys=qp->px.ys;
	int ic,ic_max=0;
	qp->intr=true;
	qp->extr=true;
	for(i=0;i<sld.nelp;i++){
		elp=sld.els[i];
		ic=poly_cross(pl, elp);
		if(ic>ic_max) ic_max=ic;

		if(ic==3) qp->bndr=true;  // all boundary
		if(ic!=1) qp->intr=false; // intersection
		if(ic!=0) qp->extr=false; // ~intersection

		if(ic >1){
			icrs=ic;
//			break;
		}
	};

	icrs=0;
	if(qp->intr) icrs=1;
	if(qp->bndr) icrs=3;

	if(lev > 6){
		(*count)++;
		//qp->icrs=ic_max;
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
			qp->chld[k]=new_QPatch(Ya,Wd);
			qp->chld[k]->par=qp;
			qp->chld[k]->lev=lev+1;
			//Qtree(qp->chld[k], cr, count);
			Qtree(qp->chld[k], sld, count);
			k++;
		}
		}
	}else{
		(*count)++;
		//qp->icrs=ic_max;
		qp->icrs=icrs;
		//qp->draw();
	};
	return(qp->lev);
};

Ellip::Ellip(){
	xc[0]=0.0; xc[1]=0.0;
	radi[0]=1.0; radi[1]=1.0;
	phi=0.0;
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
Poly::Poly(){};
Poly::Poly(int n){
	np=n;
	Poly::mem_alloc();
	xg[0]=0.0;
	xg[1]=0.0;
};
void Poly::mem_alloc(){
	xs=(double *)std::malloc(sizeof(double)*np);
	ys=(double *)std::malloc(sizeof(double)*np);
	xg[0]=0.0;
	xg[1]=0.0;
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

	return(poly_cross(Ad,Bd));
};
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

	//bbu.draw();
	//bbc.draw();
	//
	

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
