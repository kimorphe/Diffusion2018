#define DB 0
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "tcntrl.h"

/*
//---------------------------------------------------------------
//	Temperature Control
class Temp_Hist{
	public:
		double Ts,Te;	// Temperature (start,end)
		double Tdiff;	
		double Temp;	// current temperature
		bool cont_iteration; // control if another iteration will be taken
		unsigned long  nstep; 
		double p;	//exponent 0 < p < infnity 
		unsigned long istep; // current time step
		Temp_Hist(double T1, double T2, unsigned long  n_step, double p_exp);
		Temp_Hist(double T1, double T2, unsigned long n_step);
		Temp_Hist();
		double inc_Temp(); // increment temperature linearly
		double inc_Temp_exp(); // increment temperature along expolential curve
		void clear(); // initialize class object
		double tau(); 
		void load(char *fname); // load Ts, Te, nstep from a file
		int Nout;
	private:
};
*/

void Temp_Hist::load(char *fname){
	FILE *fp=fopen(fname,"r");

	char cbff[128];
	//fgets(cbff,128,fp);
	//fscanf(fp,"%d\n",&Nout);

	fgets(cbff,128,fp);
	fscanf(fp,"%le %le\n",&Ts,&Te);
	fgets(cbff,128,fp);
	fscanf(fp,"%ld\n",&nstep);
	Tdiff=Ts-Te;
	istep=0;	// current step
	cont_iteration=true;
	p=0.1;

	Temp=Ts;
	fclose(fp);

};
//---------------------------------------------------------------
Temp_Hist::Temp_Hist(){}; // default constructor
Temp_Hist::Temp_Hist(	// Constructor 1
	double T1, 	// initial temperature
	double T2, 	// final temperature
	unsigned long n_step 	// nuber of MC steps
){
	Ts=T1; 
	Te=T2;	
	Tdiff=Ts-Te;
	nstep=n_step;
	istep=0;	// current step
	cont_iteration=true;
	p=0.1;
	Temp=Ts;
};
void Temp_Hist::renew(double alph){

	Tdiff*=alph;
	Ts=Te+Tdiff;
	istep=0;
	cont_iteration=true;
	Temp=Ts;	
};

Temp_Hist::Temp_Hist(	// Constructor 2
	double T1, 	// initial temperature
	double T2, 	// final temperature
	unsigned long n_step, 	// nuber of MC steps
	double p_exp	// exponent of the rational function 
){
	Ts=T1; 
	Te=T2;	
	Tdiff=Ts-Te;
	p=p_exp;
	nstep=n_step;
	istep=0;	// current step
	cont_iteration=true;
	Temp=Ts;
};
double Temp_Hist::tau(){ // return normalized time
	// Note istep,nstep is UNSIGNED LONG, which can never be negative
	//  substitution of a negative integer results in a huger positive integer
	//return(double(istep-1)/nstep);
	return(double(istep)/nstep);
};
double Temp_Hist::inc_Temp(){	// temperature control by a rational function 
	Temp=1.0-pow((istep/(double)nstep),p);
	Temp=pow(Temp,1./p)*Tdiff+Te;
	istep++;
	if(istep >= nstep) cont_iteration=false;
	return(Temp);
};
double Temp_Hist::inc_Temp_exp(){	// exponentially decreasing temperature

	double alph=-log(Te/Ts)/nstep;
	Temp=exp(-Tdiff*istep/(double)nstep)*Ts;
	Temp=exp(-alph*istep)*Ts;
	istep++;
	if(istep >= nstep) cont_iteration=false;
	return(Temp);
};
void Temp_Hist::clear(){
	istep=0;
	cont_iteration=true;
};
//---------------------------------------------------------------
#if DB ==1
int main(){
	char f_temp[128]="temp_hist.inp";
	Temp_Hist TH;
	TH.load(f_temp);
	while(TH.cont_iteration){
		TH.inc_Temp_exp();
		printf("%ld %lf %lf\n",TH.istep,TH.tau(),TH.Temp);
		//vsl.mc(TH,false,0.6,0.1);
		//if(vsl.Etot<vsl.ncell*tol) break;
	}

	return(0);
};
#endif
