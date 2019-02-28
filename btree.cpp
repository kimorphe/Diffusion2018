#include <stdio.h>
#include <stdlib.h>

//--------------------------------------------------------
class Seg{
	public:
		double x1,x2;
		Seg(double xl, double xr);
		Seg();
		void set_xlim(double xl, double xr);
		bool is_in(double x);
		void print();
	private:
	protected:
};	
//class BinPatch: public Seg{
class BinPatch: public Seg{
	public:
		Seg A;
		BinPatch *par;
		BinPatch *chld[2];
		BinPatch();
		void print();
		void set_xlim(double xl, double xr);
		int lev;
	private:
	protected:
};
//--------------------------------------------------------
Seg::Seg(){
};
Seg::Seg(double xl, double xr){
	x1=xl;
	x2=xr;
};
void Seg::set_xlim(double xl, double xr){
	x1=xl;
	x2=xr;
};
void Seg::print(){
	//printf("(x1,x2)=(%lf, %lf)\n",x1,x2);
	printf("%lf %lf 0.0\n",x1,x2);
}
bool Seg::is_in(double x){
	bool included=true;

	if(x<x1) included=false;
	if(x>x2) included=false;

	return(included);
};
//--------------------------------------------------------
BinPatch::BinPatch(){
	par=NULL;
	chld[0]=NULL;
	chld[1]=NULL;
};
void BinPatch::print(){
//	printf("(x1,x2)=(%lf, %lf)\n",x1,x2);
	A.print();
};
void BinPatch::set_xlim(double xl,double xr){
	A.set_xlim(xl,xr);
};

BinPatch *new_Patch(double xl, double xr){
	BinPatch *bp=(BinPatch *)malloc(sizeof(BinPatch));
	bp->set_xlim(xl,xr);
	bp->chld[0]=NULL;
	bp->chld[1]=NULL;
	return(bp);
};

//--------------------------------------------------------

int set_crs(Seg A, Seg B){
	
	int icrs;
	bool a1,a2;
	a1=B.is_in(A.x1);
	a2=B.is_in(A.x2);

	if(a1 && a2 ){ // A in B
		icrs=2;
		return(icrs);
	}

	bool b1,b2;
	b1=A.is_in(B.x1);
	b2=A.is_in(B.x2);
	if(b1 && b2 ){	// B in A
		icrs=1;
		return(icrs);
	}

	if(a1||a2) return(3);

	return(0);
};

void translate_crs(int icrs){

	if(icrs==0) puts("A ^ B = phi");
	if(icrs==1) puts("B < A ");
	if(icrs==2) puts("A < B ");
	if(icrs==3) puts("A ^ B != phi");
};

int down(BinPatch *bp, Seg G,int *count){
	double xl,xr,xm;

	int lev=bp->lev;

	if(lev > 4){
		bp->print();
		(*count)++;
		return(lev);
	}

	int icrs=set_crs(G, bp->A);
	//printf("Relation of G to A is %d\n",icrs);
	//translate_crs(icrs);
	if(icrs>1){
		xl=bp->A.x1;
		xr=bp->A.x2;
		xm=(xl+xr)*0.5;
		bp->chld[0]=new_Patch(xl,xm);
		bp->chld[1]=new_Patch(xm,xr);
		bp->chld[0]->par=bp;
		bp->chld[1]->par=bp;
		bp->chld[0]->lev=lev+1;
		bp->chld[1]->lev=lev+1;
		down( bp->chld[0], G, count);
		down( bp->chld[1], G, count);
	}else{
		(*count)++;
		bp->print();
	};
	return(bp->lev);
};
//--------------------------------------------------------
void leaves(BinPatch *bp,int *count,BinPatch *bps){
	if(bp->chld[0]==NULL){
		bp->print();
		printf("idx=%d\n",*count);
		bps[*count]=(*bp);
		(*count)++;
	}else{
		leaves(bp->chld[0],count,bps);
		leaves(bp->chld[1],count,bps);
	};
};
//--------------------------------------------------------
int main(){
	int count;

	Seg G(0.3,0.7);

	BinPatch bp0;	// root node
	bp0.A.set_xlim(0.0,1.0);
	//bp.print();
	bp0.lev=0;

	printf("# G=");
	G.print();
	count=0;
	down(&bp0,G,&count);
	printf("count=%d\n",count);

	puts("-----------------------------------");
	//count=0;

	BinPatch *bps=(BinPatch *)malloc(sizeof(BinPatch)*count);
	int idx;
	leaves(&bp0,&idx,bps);

	puts("-----------------------------------");
	for(idx=0;idx<count;idx++){
		printf("Lev=%d, ",bps[idx].lev);
		bps[idx].print();
	};


	return(0);
};
