#include<stdio.h>
#include<stdlib.h>

int main(int argc, char *argv[]){
	int ndir=10;
	char  cbff[128],fname[128];
	FILE *fp;
	double *u2,*ux,*uy;
	double v1,v2,v3;

	int i,j,k,ndat;

	for(k=1;k<=8;k++){
		printf("Sr%d\n",k);
	for(i=1;i<=ndir;i++){
		sprintf(fname,"Sr%d/SEED%d/u2b.out",k,i);
		printf("%s\n",fname);
		if(i==1){
			ndat=0;
			fp=fopen(fname,"r");
				while(fgets(cbff,128,fp)!=NULL) ndat++;	
			fclose(fp);
			u2=(double *)malloc(sizeof(double)*ndat);
			ux=(double *)malloc(sizeof(double)*ndat);
			uy=(double *)malloc(sizeof(double)*ndat);

			printf("ndat=%d\n",ndat);
			for(j=0;j<ndat;j++){
				u2[j]=0.0;
				ux[j]=0.0;
				uy[j]=0.0;
			}
		}

		fp=fopen(fname,"r");
		for(j=0;j<ndat;j++){
			fscanf(fp,"%lf %lf %lf\n",&v1,&v2,&v3);
			u2[j]+=v1;
			ux[j]+=v2;
			uy[j]+=v3;
		}
		fclose(fp);
	}	

	sprintf(fname,"Sr%d/u2b.out",k);
	fp=fopen(fname,"w");
	for(j=0;j<ndat;j++){
		fscanf(fp,"%lf %lf %lf\n",&v1,&v2,&v3);
		u2[j]/=ndir;
		ux[j]/=ndir;
		uy[j]/=ndir;
		fprintf(fp,"%lf %lf %lf\n",u2[j],ux[j],uy[j]);
	}
	fclose(fp);

	free(u2);
	free(ux);
	free(uy);
	}

	return(0);
};
