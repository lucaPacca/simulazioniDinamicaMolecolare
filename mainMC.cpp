#include <iostream>
#include <sstream> 
#include "SistemaMRT.h"
using namespace std;

int main()
{
	//Nomi dei file su cui salvare i dati
	string path = "./";
	string dati_eq = path+"Dati_eq_gas.csv";
	string dati_stat = path+"Dati_stat_gas.csv";
	string dati_oss = path + "Dati_oss_gas.csv";
	string dati_cv = path + "Dati_cv_gas.csv";
	string dati_err = path + "Dati_err_gas.csv";
	double T= 1.1; //Temperatura adimensionale
	double rho = 0.01;//Densità adimensionale
	double Delta = 0.01;
	double costanti_d[3] = {T,rho,Delta};//Vettore costanti double
    int n= 8; //Numero particelle per scatola 
	int N_eq = 50000; //Numero di passi per la fase di equilibriatura
	int N_stat = 1000000; //Numero di passi per la fase statistica
	int N_cv = 2000;
	int M = 1; //Numero di particelle per cella
	int costanti_i[5] = {n,N_eq,N_stat,N_cv,M}; //Vettore costanti int
	SistemaMRT s(costanti_d, costanti_i);
	s.CreaSistema();
	s.print_condizioni_iniziali();
	s.fase_equilibriatura(&dati_eq);
	s.fase_statistica(&dati_stat);
	s.calcola_osservabili(&dati_oss,&dati_cv);
	s.block(&dati_err);
	return 0;	
}



