#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <random>
#include "set2d.h"
#include "domain.h"

int main(){
	char fname[128]="kcell.dat";
	Circ cdat;

	int nx=200;
	int ny=100;
	Dom2D dom(nx,ny);

	dom.Xa[0]=0.0; dom.Xa[1]=0.0;
	dom.Xb[0]=2.0; dom.Xb[1]=1.0;
	dom.set_dx();
		cdat.xc[0]=1.0;
		cdat.xc[1]=0.0;
		cdat.radi=0.5;
	dom.perfo0(cdat,1);
		cdat.xc[0]=1.0;
		cdat.xc[1]=1.0;
		cdat.radi=0.5;
	dom.perfo0(cdat,1);
		cdat.xc[0]=2.0;
		cdat.xc[1]=1.0;
		cdat.radi=.25;
	dom.perfo0(cdat,1);

	dom.out_kcell(fname);

	return(0);
};
