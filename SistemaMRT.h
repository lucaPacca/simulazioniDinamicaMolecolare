#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdlib.h>
#include <fstream>

using namespace std;
//Costanti
#define pi  3.141592653589793

/*
	Libreria che svolge simulazioni e misura di osservabili
	dati parametri con il metodo di Monte Carlo.
	
	Classi e funzioni:
	
	LJ							 44 - POTENZIALE LENNARD JONES
	DLJ							 53 - DERIVATA POTENZIALE
	IntegraleSimpson			 62 - INTEGRAZIONE SIMPSON PER UNA MATRICE
 
 Posizioni						 78 - CLASSE POSIZIONI
  public
  
 SistemaMRT						 92 - CLASSE DEL SISTEMA 
  public:
	CreaSistema					135 - INIZIALIZZA LE VARIABILI DEL SISTEMA
	print_condizioni_iniziali	146 - STAMPA LE OSSERVABILI A t=0
	fase_equilibriatura			156 - FASE DI EQUILIBRIATURA, LA DURATA DIPENDE DA N_EQ
	fase_statistica				214 - FASE DI STATISTICA, LA DURATA DIPENDE DA N_STAT
	calcola_osservabili			259 - CALCOLA LE OSSERVABILI DEL SISTEMA DOPO LA FASE STATISTICA
	block						306	-  CALCOLA CORRELAZIONI SU BLOCCHI DI ITERAZIONI
  
  private:
	cpc							394	- TAGLIA DISTANZE SECONDO LA COND. PER. CONTORNO
	init_r						401 - CREA IL RETICOLO
	calcE						448 - CALCOLA ENERGIA E LAVORO TRA PARTICELLE
	passo_metropolis			474 - AGGIORNA POSIZIONI CON M(RT)^2
	update_cv					504	- AGGIORNA LA FUNZIONE DI CORRELAZIONE
	norm_cv						517	- NORMALIZZA LA FUNZIONE DI CORRELAZIONE
	calc_tau					541 - CALCOLA GLI LA LUNGHEZZA DI CORRELAZIONE
*/

//----------------- POTENZIALE LENNARD JONES
double LJ(double x2) 
{
	double v1 = x2 * x2 * x2;
	double v2 = v1 * v1;
	return 4 * (1 / v2 - 1 / v1);
}

//----------------- DERIVATA POTENZIALE

double DLJ(double x)
{
	double v1 = x * x * x * x * x * x;
	double v2 = v1 * v1;
	return 24 * (-2 / (v2 * x) + 1 / (v1 * x));
}

//----------------- INTEGRAZIONE SIMPSON PER UNA MATRICE
 
double IntegraleSimpson(double* f, double h, int N) {
	double result = (*f) + (*(f + N - 1));
	for (int i = 1; i < N - 1; i++) {
		f++;
		if (i % 2 == 0) {
			result += (*f) * 2;
		}
		else {
			result += (*f) * 4;
		}
	}
	return result * h / 3;
}

//----------------- CLASSE POSIZIONI

class Posizioni {
    public:
		double** s;
		Posizioni(int N)
		{
			s = new double* [N];
			for (int i = 0; i < N; i++)
			{
				s[i] = new double[3];
			}
	    }
};

//----------------- CLASSE DEL SISTEMA 
class SistemaMRT {
public:

	//----------------- CLASSE DEL SISTEMA 

	SistemaMRT(double* costanti_d, int* costanti_i)
	{
		//Dichiaro le costanti
		this->T = *(costanti_d);
		this->rho = *(costanti_d + 1);
		this->Delta = *(costanti_d + 2);

		this->n = *costanti_i;
		this->N_eq = *(costanti_i + 1);
		this->N_stat = *(costanti_i + 2);
		this->N_cv = *(costanti_i + 3);
		this->M=*(costanti_i + 4);

		N = n * n * n * M; //Numero di particelle
		L = pow(N / rho, 1. / 3); //Lato della scatola
		V = L * L * L;
		dDelta = Delta * 0.5;
		EM = new double[N_cv];
		E2M = new double[N_cv];
		WM = new double[N_cv];
		CE = new double[N_cv];
		CE2 = new double[N_cv];
		CW = new double[N_cv];
		for (int i = 0; i < N_cv; i++)
		{
			CE[i] = 0, CE2[i] = 0, CW[i] = 0;
			
		}
		s1 = new Posizioni(N);
		s2 = new Posizioni(N);
		temp = new Posizioni(N);
		Wblock = new double[N_stat];
		Eblock = new double[N_stat];

	}
	
	//----------------- INIZIALIZZA LE VARIABILI DEL SISTEMA
	
	void CreaSistema()
	{
		init_r();
		calcE(s1);
		E = E_new;
		W = W_new;
		test[0] = 2;
	}
	
	//----------------- STAMPA LE OSSERVABILI A t=0
	
	void print_condizioni_iniziali() 
	{
		cout << " DINAMICA MOLECOLARE : Simulazione Monte Carlo NVT \n" << endl;
		cout << " Temperatura: " << T << endl;
		cout << " Energia potenziale iniziale: " << E << endl;
		cout << " Viriale iniziale: " << W << endl;
	}
	
	//----------------- FASE DI EQUILIBRIATURA, LA DURATA DIPENDE DA N_EQ
	
	void fase_equilibriatura(string* filenameeq = NULL)
	{
		ofstream fileeq;
		if (filenameeq != NULL)
		{
			fileeq.open(*filenameeq);
		}
		cout << "\nInizio fase di equilibriatura\n" << endl;
		for (int i = 0; i < N_eq; i++)
		{
			if (passo_metropolis())
			{
				acc++;
				acc_tot++;
			}

			if (i % Step_control == 0 && i< N_eq/2&&i != 0) {
				fraz_acc = (double)acc / Step_control;

				if (fraz_acc > acc_max) {
					Delta += dDelta;
					test[1] = 1;
				}
				else if (fraz_acc < acc_min) {
					Delta -= dDelta;
					test[1] = 0;
				}

				if (Delta <= 0) {
					Delta += dDelta;
					dDelta *= 0.5;

				}
				else if (test[0] + test[1] == 1)
				{
					dDelta *= 0.5;
				}
				test[0] = test[1];
				acc = 0;
			}
			
			
			if (i % Step_print == 0) {
				cout << "i: " << i << " U: " << E <<" Viriale: "<<W<< " Delta: " << Delta << " Acc: " << (double)acc_tot / (double)(i+1) << endl;
			}
			if (i % Step_save == 0 && filenameeq != NULL)
			{
				fileeq << i << ";" << E << ";" <<W<<";" << Delta << ";" << (double)acc_tot/(double)(i+1) << endl;
			}
		}
		fileeq.close();
		cout << "\nFine fase di equilibriatura\n" << endl;
		cout << "Acceptance rate :" << (double)acc_tot / (double)N_eq << endl;
		cout << "Delta :" << Delta << endl;
	}
	
	//----------------- FASE DI STATISTICA, LA DURATA DIPENDE DA N_STAT
	
	void fase_statistica(string* filenamestat = NULL)
	{
		ofstream filestat;
		if (filenamestat != NULL)
		{
			filestat.open(*filenamestat);
		}
		cout << "\nInizio fase di statistica\n" << endl;
		acc_tot = 0,idx=0;
		for (int i = 0; i < N_stat; i++)
		{
			if (passo_metropolis()) {
				
				 acc_tot++;
			}
			idx = (idx + 1) % N_cv;
			EM[idx] = E;
			WM[idx] = W;
			
			
			Eblock[i] = E;
			Wblock[i] = W;
			if (i >= N_cv)
			{
				update_cv();
			}
			Esum += E;
			E2sum += E2;
			Wsum += W;
			W2sum += W * W;
			 if (i % Step_print == 0) {
				 cout << "i: " << i << " U: " << E << " W: " << W << " Acceptance rate :" << (double)acc_tot / (double)(i+1) << endl;
			 }
			 if (i % Step_save == 0 && filenamestat != NULL)
			 {
				 filestat << i << ";" << E << ";" << W << ";"  << (double)acc_tot / (double)(i + 1) << endl;
			 }
		}
		filestat.close();
		cout << "\nFine fase di statistica\n" << endl;
		cout << "Acceptance rate :" << (double)acc_tot / (double)N_stat << endl;
	}
	
	//----------------- CALCOLA LE OSSERVABILI DEL SISTEMA DOPO LA FASE STATISTICA
	
	void calcola_osservabili(string* filenameoss = NULL, string* filenamecv = NULL) {
		
		double Emean = Esum / N_stat;
		double E2mean = E2sum / N_stat;
		double Wmean = Wsum / N_stat;
		double W2mean = W2sum / N_stat;
		double Evar = E2sum / N_stat - pow(Emean, 2);
		double Estd = sqrt(Evar / N_stat);
		double Wvar = W2sum / N_stat - pow(Wmean, 2);
		double Wstd = sqrt(Wvar / N_stat);

		double pres_mean = Wmean / (3. * V) + rho * T;
		double pres_std = Wstd / (3. * V);
	
		norm_cv(Emean,Wmean,N_cv, N_stat - N_cv + 1,filenamecv);
		double tau_E, tau_W;
		tau_E=calc_tau(CE,N_cv);
		tau_W=calc_tau(CW,N_cv);
		
		double Estdcorr = sqrt(((E2sum / (double)N_stat - Emean * Emean) /((double)N_stat)) * (tau_E+1));
		double Wstdcorr = sqrt(((W2sum / (double)N_stat - Wmean * Wmean) / ((double)N_stat)) * (tau_W+1));
		double pres_stdcorr = Wstdcorr / (3. * V);
		

		cout << "\nEnergia potenziale: " << Emean << " +- " << Estd << " - " << Estdcorr <<" tau: "<<tau_E<< endl;
		cout << "Viriale: " << Wmean << " +- " << Wstd << " - " << Wstdcorr << " tau: " << tau_W<< endl;
		cout << "Pressione: " << pres_mean << " +- " << pres_std << " - " << pres_stdcorr<< endl;
		
		
		ofstream fileoss;
		if (filenameoss != NULL)
		{
			fileoss.open(*filenameoss);
		}
		fileoss << rho << endl;
		fileoss << Emean << ";" << Estd << ";" << Estdcorr << endl;
		fileoss << Wmean << ";" << Wstd << ";" << Wstdcorr << endl;
		fileoss << pres_mean << ";" << pres_std << ";" << pres_stdcorr << endl;
		fileoss << tau_E << ";" << tau_W << endl;
		fileoss << Delta << ";" << (double)acc_tot / (double)N_stat << endl;

		fileoss.close();

	}
	
	//----------------- CALCOLA CORRELAZIONI SU BLOCCHI DI ITERAZIONI
	
	void block(string* filenameerr=NULL)
	{
		int Bck[110];
		int value = 2;
		for (int i = 0; i < 110; i++)
		{
			Bck[i] = value;
			if (value == 1024)
			{
				value = 2000;
			}
			else if (value < 1024)
			{
				value *= 2;
			}
			else
			{
				value += 2000;
			}
			cout << Bck[i] << endl;
		}
		ofstream fileerr;
		if (filenameerr != NULL)
		{
			fileerr.open(*filenameerr);
		}
		double* Wsumblock,*Esumblock;
		double Waverage, Eaverage;
		double Wstd , Estd;
		double Werr[110], Eerr[110];
		int index;
		for (int i = 0; i < 110; i++)
		{
			Wstd = 0, Estd = 0, Waverage = 0, Eaverage = 0;
			index = 0;
			Wsumblock = new double[(int)(N_stat / Bck[i])];
			Esumblock = new double[(int)(N_stat / Bck[i])];
			for (int j = 0; j < (int)(N_stat/Bck[i]); j++)
			{
				for (int k = 0; k < Bck[i]; k++)
				{
					Wsumblock[j] += Wblock[index];
					Esumblock[j] += Eblock[index];
					index++;
				}
				Wsumblock[j] /= Bck[i];
				Esumblock[j] /= Bck[i];
				Waverage += Wsumblock[j];
				Eaverage += Esumblock[j];

			}
			Waverage /= (N_stat / Bck[i]);
			Eaverage /= (N_stat / Bck[i]);
			for (int j = 0; j < (int)(N_stat/Bck[i]); j++)
			{
				Wstd += (Wsumblock[j] - Waverage) * (Wsumblock[j] - Waverage);
				Estd += (Esumblock[j] - Eaverage) * (Esumblock[j] - Eaverage);
			}
			Wstd /= (N_stat / Bck[i]);
			Estd /= (N_stat / Bck[i]);
			Werr[i] = sqrt(Wstd / (N_stat / Bck[i]));
			Eerr[i] = sqrt(Estd / (N_stat / Bck[i]));
			free(Wsumblock);
			free(Esumblock);	
		}
		for (int i = 0; i < 110; i++)
		{
			fileerr<< Eerr[i] <<";" << Werr[i] <<";"<< Bck[i] << endl;
		}
		fileerr.close();
	}
	
//----------------- PARTE PRIVATA
	
private:
	int N,n,N_eq,N_stat,N_cv,M;
	int acc=0,acc_tot=0,idx;
	int Step_print = 1000, Step_control = 200, Step_save = 10, test[2];
	double rho, T, L, V, E, E2,W , Esum = 0, E2sum = 0, E4sum = 0, Wsum = 0, W2sum = 0, dDelta,Delta, fraz_acc;
	double acc_max = 0.55, acc_min = 0.45,E_new,W_new;
	double* EM, * E2M, * WM, * CE, * CE2, * CW;
	Posizioni* s1;
	Posizioni* s2;
	Posizioni* temp;
	double* Wblock, * Eblock;
	
	//----------------- CONDIZIONE PERIODICA AL CONTORNO
	
	double cpc(double s) 
	{
		return s - rint(s);
	}
	
	//----------------- CREA IL RETICOLO
	
	void init_r()
	{
		double dL = 1./n;
		int i = 0;
		double BCC[2][3] = { {0,0,0},{0.5,0.5,0.5} };
		double FCC[4][3] = { {0,0,0},{0.5,0.5,0},{0.5,0,0.5},{0,0.5,0.5} };
		for (int x = 0; x < n; x++)
		{
			for (int y = 0; y < n; y++)
			{
				for (int z = 0; z < n; z++)
				{
					if (M == 1)
					{
						s1->s[i][0] = -0.5 + dL * x;
						s1->s[i][1] = -0.5 + dL * y;
						s1->s[i][2] = -0.5 + dL * z;
						i++;
					}
					if (M == 2)
					{
						for (int j = 0; j < M; j++)
						{
							s1->s[i][0] = -0.5 + dL * ((double)x + BCC[j][0]);
							s1->s[i][1] = -0.5 + dL * ((double)y + BCC[j][1]);
							s1->s[i][2] = -0.5 + dL * ((double)z + BCC[j][2]);
							i++;

						}
					}
					if (M == 4)
					{
						for (int j = 0; j < M; j++)
						{
							s1->s[i][0] = -0.5 + dL * ((double)x + FCC[j][0]);
							s1->s[i][1] = -0.5 + dL * ((double)y + FCC[j][1]);
							s1->s[i][2] = -0.5 + dL * ((double)z + FCC[j][2]);
							i++;
						}
					}
				}
			}
		}
	}
	
	//----------------- CALCOLA ENERGIA E LAVORO TRA PARTICELLE
	
	void calcE(Posizioni *p)
	{
		double ds2,ds_mod;
		E_new = 0, W_new = 0;
		for (int i = 1; i < N; i++)
		{
			for (int j = 0; j < i; j++)
			{
				ds2 = 0;
				for (int k = 0; k < 3; k++)
				{
					ds2 += pow(cpc(p->s[i][k] - p->s[j][k]), 2);
				}
				ds_mod = sqrt(ds2);
				if (ds_mod < 0.5) {

					E_new += LJ(ds2 * L * L);
					W_new -= DLJ(ds_mod * L) * ds_mod*L;
				}

			}
		}
	}
	
	//----------------- AGGIORNA POSIZIONI CON M(RT)^2
	
	bool passo_metropolis()
	{
		double si, prob;
		
		for (int i = 0; i < N; i++)
		{
			for (int k = 0; k < 3; k++)
			{
				si = (double)rand() / (RAND_MAX + 1.);
				s2->s[i][k] = cpc(s1->s[i][k] + Delta * (si - 0.5));
			}
		}
		calcE(s2);
		prob = exp(-(E_new - E) / T);
		si = (double)rand() / (RAND_MAX + 1.);
		if (si < prob)
		{
			temp = s1;
			s1 = s2;
			s2 = temp;
			E = E_new;
			E2 = E*E;
			W = W_new;
			return true;
		}
		return false;
	}
	
	//----------------- AGGIORNA LA FUNZIONE DI CORRELAZIONE
	
	void update_cv()
	{
		int idx2;
		for (int i = 0; i <N_cv; i++) {
			idx2 = (idx - i + N_cv) % N_cv;
			CE[i] += EM[idx] * EM[idx2];
			//CE2[i] += E2M[idx] * E2M[idx2];
			CW[i] += WM[idx] * WM[idx2];
		}
	}
	
	//----------------- NORMALIZZA LA FUNZIONE DI CORRELAZIONE
	
	void norm_cv(double meanE, double meanW, int len, int N_step,string* filenamecv=NULL)
	{
		ofstream filecv;
		if (filenamecv != NULL)
		{
			filecv.open(*filenamecv);
		}
		double meanE2 = meanE * meanE;
		double meanW2 = meanW * meanW;
		CE[0] = CE[0] / (double)N_step - meanE2;
		CW[0] = CW[0] / (double)N_step - meanW2;
		for (int i = 1; i < len; i++) {
			CE[i] = ((CE[i] / (double)N_step) - meanE2) / CE[0];
			CW[i] = ((CW[i] / (double)N_step) - meanW2) / CW[0];
			filecv << CE[i] << ";" << CW[i] << ";" << i << endl;
		}
		CE[0] = 1;
		CW[0] = 1;
		filecv << CE[0] << ";" << CW[0] << ";" << 0 << endl;
		filecv.close();
	}
	
	//----------------- CALCOLA GLI LA LUNGHEZZA DI CORRELAZIONE
	
	double calc_tau(double* C,int len)
	{
		double tau = 0, inv_e = exp(-1);
		double c = 1;
		while (tau == 0)
		{
			for (int i = 1; i < len; i++) {
				if (C[i] <= inv_e && tau == 0) {

					if (i != 1)
					{
						tau = (inv_e - C[i]) / (C[i] - C[i - 1]) + (double)i*c;
					}
					else
					{
						tau = (inv_e - C[i]) / (C[i] - 1) + (double)i*c;
					}
				}
			}
			if (tau == 0)
			{
				c++;
				inv_e = exp(-1. / c);
			}
			if (c == 10)
			{
				cout << "Errore,la correlazione non converge" << endl;
				tau = -1;
			}
		}
		return tau;
	}
	
};
