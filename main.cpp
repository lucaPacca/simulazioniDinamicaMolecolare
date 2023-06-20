#include <iostream>
#include "SistemaLJ.h"
using namespace std;

int main()
{
	//Nomi dei file su cui salvare i dati
	string path = "C:/Users/Federico/Desktop/Fisica computazionale/Data/";
	string dati_eq_gas = path+"Dati_eq_liq.csv";
	string dati_stat_gas = path+"Dati_stat_liq.csv";
	string dati_oss_gas = path + "Dati_oss_liq.csv";
	string dati_gr_gas = path + "Dati_gr_liq.csv";
	string dati_cv_gas = path + "Dati_cv_liq.csv";

	double T= 1.1; //Temperatura adimensionale
	double rho = 0.01;//Densità adimensionale
	double dt = 0.001; //dt adimensionale
	double costanti_d[3] = {T,rho,dt};//Vettore costanti double
    int n = 5; //Numero particelle per lato 
	int N_eq = 50000; //Numero di passi per la fase di equilibriatura
	int N_stat = 50000; //Numero di passi per la fase statistica
	int N_bin = 10000; //Numero di bin per il calcolo di g(r)
	int M = 1; //Numero di particelle per cella(1 per CC,2 per BCC,4 per FCC)
	int N_cv = 5000; //Lunghezza del vettore cv
	int costanti_i[6] = {n,N_eq,N_stat,N_bin,M,N_cv}; //Vettore costanti int

	SistemaLJ s(costanti_d, costanti_i);
	s.CreaSistema();
	s.print_condizioni_iniziali();
	s.fase_equilibriatura(&dati_eq_gas); //Se voglio salvare su file: s.fase_equilibriatura(&dati_eq_gas);
	s.fase_statistica(&dati_stat_gas);
	s.calcola_osservabili(&dati_oss_gas);
	s.save_Gr(&dati_gr_gas);
	s.save_cv(&dati_cv_gas);

	/*
	Oppure:

	SistemaLJ* s = new SistemaLJ(costanti_d, costanti_i);
	s->CreaSistema();
	s->print_condizioni_iniziali();
	s->fase_equilibriatura();
	s->fase_statistica();
	s->calcola_osservabili();
	s->save_cv();
	delete s;
	*/
	return 0;
}


