#include<stdio.h>
#include<stdlib.h>
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
		void renew(double alph);
	private:
};

