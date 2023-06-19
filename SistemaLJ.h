﻿#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdlib.h>
#include <fstream>

using namespace std;

//Costanti
#define pi  3.141592653589793


/*
	Libreria che svolge simulazioni e misura di osservabili
	dati parametri.
	
	Classe (SistemaLJ) e funzioni:
	
	LJ							 47 - POTENZIALE LENNARD JONES
	DLJ							 55 - DERIVATA POTENZIALE
	IntegraleSimpson			 64 - INTEGRAZIONE SIMPSON PER UNA MATRICE
	
  public:
	CreaSistema					120 - INIZIALIZZA LE VARIABILI DEL SISTEMA
	print_condizioni_iniziali	131 - STAMPA LE OSSERVABILI A t=0
	fase_equilibriatura			142 - FASE DI EQUILIBRIATURA, LA DURATA DIPENDE DA N_EQ
	fase_statistica				169 - FASE DI STATISTICA, LA DURATA DIPENDE DA N_STAT
	save_cv						225 - SALVA LA FUNZIONE DI AUTOCORRELAZIONE
	save_Gr						241 - SALVA LA FUNZIONE DISTRIBUZIONE DI COPPIA
	calcola_osservabili			254 - CALCOLA LE OSSERVABILI DEL SISTEMA DOPO LA FASE STATISTICA
  
  private:
	init_cv						322 - INIZIALIZZA A 0 LA FUNZIONE DI CORRELAZIONE
	update_cv					332 - AGGIORNA LA FUNZIONE DI CORRELAZIONE
	init_r						365 - CREA IL RETICOLO
	init_v						453 - CREA LA DISTRIBUZIONE DELLE VELOCITÀ
	cpc							474 - CONDIZIONE PERIODICA AL CONTORNO
	calcF						481 - CALCOLA LE FORZE AGENTI SU OGNI PARTICELLA
	passo_verlet				525 - AGGIORNA LE POSIZIONI DELLE PARTICELLE 
	init_Gr						556 - CREA IL VETTORE DI DISTRIBUZIONE DI COPPIA
	
*/


//----------------- POTENZIALE LENNARD JONES

double LJ(double x2) 
{
	double v1 = x2 * x2 * x2;
	double v2 = v1 * v1;
	return 4 * (1 / v2 - 1 / v1);

//----------------- DERIVATA POTENZIALE

double DLJ(double x) 
{
	double v1 = x * x * x * x * x * x;
	double v2 = v1 * v1;
	return 24 * (-2 / (v2 * x) + 1 / (v1 * x));
}

//-----------------INTEGRAZIONE SIMPSON PER UNA MATRICE

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

//----------------- CLASSE DEL SISTEMA 

class SistemaLJ {
public:
	SistemaLJ(double* costanti_d, int* costanti_i)
	{

		//Dichiaro le costanti
		
		this->T = *(costanti_d );
		this->rho = *(costanti_d + 1);
		this->dt = *(costanti_d + 2);

		this->n = *costanti_i;					//Numero particelle per lato 
		this->N_eq = *(costanti_i + 1);			//Numero di passi per la fase di equilibriatura
		this->N_stat = *(costanti_i + 2);		//Numero di passi per la fase statistica
		this->N_bin = *(costanti_i + 3);		//Numero di bin per il calcolo di g(r)
		this->M=*(costanti_i + 4);				//Numero di particelle per cella
		this->N_cv = *(costanti_i + 5);			//Lunghezza del vettore cv
		N = n * n * n * M; 						//Numero di particelle
		L = pow(N/rho, 1. / 3); 				//Lato della scatola
		L2 = L / 2;								//Mezzo lato cella
		dR = L2 / N_bin;  						//Intervallo per il calcolo di g(r)
		V = L * L * L;							//Volume 
		gr = new double[N_bin];					//Contatore distanze particelle
		norm = new double[N_bin];				//valore distanze particelle
		Cv = new double[N_cv];					//volume intervallo distanza tra particelle
		r = new double* [N];					//Insieme posizioni
		F = new double* [N];					//Insieme forze
		v = new double* [N];					//Insieme velocità
		v0 = new double* [N];
		for (int i = 0; i < N; i++)
		{
			r[i] = new double[3];
			v[i] = new double[3];
			v0[i] = new double[3];
			F[i] = new double[3];
		}
	}

	//----------------- INIZIALIZZA LE VARIABILI DEL SISTEMA

	void CreaSistema() 
	{
		t = 0;
		init_v();
		init_r();
		init_cv();
		init_Gr();
	}

	//----------------- STAMPA LE OSSERVABILI A t=0

	void print_condizioni_iniziali() 
	{
		cout << " DINAMICA MOLECOLARE : fluido di Lennard - Jones \n" << endl;
		cout << " Energia cinetica iniziale: " << E[0] << endl;
		cout << " Energia potenziale iniziale: " << E[1] << endl;
		cout << " Energia totale iniziale: " << E[2] << endl;
		cout << " Temperatura iniziale: " << Ts << endl;
	}

	//----------------- FASE DI EQUILIBRIATURA, LA DURATA DIPENDE DA N_EQ

	void fase_equilibriatura(string* filename = NULL) 
	{
		ofstream file;
		if (filename != NULL)
		{
			file.open(*filename);
		}
		cout << "Inizio fase di equilibriatura\n" << endl;
		for (int i = 0; i < N_eq; i++)
		{
			if (i % Step_print == 0)
			{
				cout << "t: " << t << " K: " << E[0] << " U: " << E[1] << " E: " << E[2] << " T: " << Ts << " ";
				cout << (double)i * 100 / N_eq << "%" << endl;
			}
			if (i % Step_save == 0 && Step_save != 0 && filename != NULL)
			{
				file <<setprecision(8)<<t << ";" << E[0] << ";" << E[1] << ";" << E[2] << endl;
			}
			passo_verlet();
		}
		file.close();
		cout << "\nFine fase di equilibriatura\n" << endl;
	}

	//----------------- FASE DI STATISTICA, LA DURATA DIPENDE DA N_STAT

	void fase_statistica(string* filename2 = NULL) 
	{
		cout << "Inizio fase statistica\n" << endl;
		fase_stat = true;
		ofstream file;

		if (filename2 != NULL)
		{
			file.open(*filename2);
		}
		
		for (int i = 0; i < N_stat; i++)
		{
			
			Tstat++;
			if (i % Step_print == 0)
			{
				cout << "t: " << t << " K:" << E[0] << " U:" << E[1] << " E:" << E[2] << " T:" << Ts << " ";
				cout <<(double)i * 100 / N_stat << "%" << endl;
			}
			if (i % Step_save == 0 && filename2 != NULL)
			{
				file << setprecision(8) << t << ";" << E[0] << ";" << E[1] << ";" << E[2] << ";" << Ts << ";" << pres_i << endl;
			}
			passo_verlet();
			update_cv();
		
			for (int j = 0; j < 3; j++)
			{
				E_sum[j] += E[j];
				E2_sum[j] += E[j]*E[j];
				
			}
			viriale_sum += viriale;
			viriale_sum2 += pow(viriale, 2);
			
		}
		for (int i = 0; i < N_bin; i++)
		{
			gr[i] /= norm[i] * N_stat;
			
		}
		for (int i = 1; i < N_cv; i++)
		{
			Cv[i] = Cv[i] / Cv[0];
		}
		if (filename2 != NULL)
		{
			file.close();
		}
		filetr.close();
	
	}

	//----------------- SALVA LA FUNZIONE DI AUTOCORRELAZIONE

	void save_cv(string* filenamecv) 
	{
		double D;
		
		ofstream filecv(*filenamecv);
        for (int i = 2; i <= N_cv; i++)
		{
		    D = 2*IntegraleSimpson(Cv, dt, i);
			filecv << setprecision(8) << Cv[i] << ";" << D << ";" << i * dt << endl;
            
		}
		filecv.close();
	}

	//----------------- SALVA LA FUNZIONE DISTRIBUZIONE DI COPPIA

    void save_Gr(string* filenamegr) 
	{
		ofstream filegr(*filenamegr);
		for (int i = 0; i < N_bin; i++)
		{
			dist = (i + 0.5) * dR;
			filegr << setprecision(8) << dist << ";" << gr[i] << endl;
		}
		filegr.close();
	}

	//----------------- CALCOLA LE OSSERVABILI DEL SISTEMA DOPO LA FASE STATISTICA

	void calcola_osservabili(string* filenameoss = NULL) 
	{
		for (int j = 0; j < 3; j++)
		{
			E_mean[j] = E_sum[j] / N_stat;
			E_var[j] = E2_sum[j] / N_stat - pow(E_mean[j], 2);
			E_std[j] = sqrt(E_var[j] / N_stat);
			
		}
		double T_mean= 2 * E_mean[0] / (3 * N);
		double viriale_mean = viriale_sum / N_stat;
		double viriale_std= sqrt((viriale_sum2 / N_stat - pow(viriale_mean, 2))/N_stat);
		double pres_mean = rho * T_mean + viriale_mean / (3 * V);
		double pres_std = sqrt((viriale_sum2 / N_stat - pow(viriale_mean, 2)) / N_stat) / (3 * V);

		cout << "\nEnergia cinetica: " << E_mean[0] << " +- " << E_std[0] << endl;
		cout << "Energia potenziale: " << E_mean[1] << " +- " << E_std[1] << endl;
		cout << "Energia totale: "<<setprecision(10) << E_mean[2] << " +- " << E_std[2] << endl;
		cout << "Viriale: " << viriale_mean << " +- " << viriale_std << endl;
		cout << "Pressione: "  << pres_mean << " +- " << pres_std << endl;
		cout << "Temperatura: "  << 2*E_mean[0]/(3*N)<< " +- " <<2*E_std[0]/(3*N) << endl;
		double* INT = new double[N_bin];
		for (int i = 0; i < N_bin; i++)
		{
			dist = (i + 0.5) * dR;
			INT[i] = gr[i] * dist * dist * LJ(dist*dist);
		}
		double dV2 = IntegraleSimpson(INT, dR, N_bin) * rho * ((N - 1)) * pi*2;
		double dV = -(3 * pow(L2, 6) + 1) / (9 * pow(L2, 9)) * rho * ((N - 1)) * 4 * pi;
		double dW = (30 * pow(L2, 6) - 24) / (5 * pow(L2, 10)) * rho * ((N - 1)) * 4 * pi;
		double pres_cor = rho * T + (viriale_mean+dW) / (3 * V);
		cout << "Correzione energia: " << dV << endl;
		cout << "Correzione energia(g(r)): " << dV2 << endl;
		cout << "Energia potenziale corretta: " << E_mean[1]+dV << endl;
		cout << "dW: " << dW << endl;
		cout << "dP: " << pres_cor << endl;
		
		ofstream fileoss;
		if (filenameoss != NULL)
		{
			fileoss.open(*filenameoss);
		}
		fileoss << setprecision(8) << E_mean[0] << ";" << E_std[0] << endl;
		fileoss << setprecision(8) << E_mean[1] << ";" << E_std[1] << endl;
		fileoss << setprecision(8) << E_mean[2] << ";" << E_std[2] << endl;
		fileoss << setprecision(8) << viriale_mean << ";" << viriale_std << endl;
		fileoss << setprecision(8) << pres_mean << ";" << pres_std << endl;
		fileoss << setprecision(8) << 2 * E_mean[0] / (3 * N) << ";" << 2 * E_std[0] / (3 * N) << endl;
		fileoss << setprecision(8) << dV << ";" << dV2 << ";" << dW <<";" << pres_cor << endl;
		if (filenameoss != NULL)
		{
			fileoss.close();
		}
	}

//----------------- PARTE PRIVATA DELLA CLASSE

private:
	int N,N_eq ,N_stat,M,Step_print = 1000, Step_save = 10,N_bin,N_cv;
	double t, E[3], E_mean[3], E_var[3], E_std[3], E_sum[3]={0,0,0}, E2_sum[3]={0,0,0}, viriale, n, rho, T, dt, L, L2, Ts,V,T_sum = 0, T_sum2 = 0, viriale_sum = 0, viriale_sum2 = 0, pres_mean, pres_std, pres_i, dR, dist;
	double* gr, * norm,*Cv;
	double v0cm[3],vcm[3];
	int Tstat = -1;
	bool fase_stat = false;
	double** r, ** F,**v,**v0;

	//----------------- INIZIALIZZA A 0 LA FUNZIONE DI CORRELAZIONE

	void init_cv()
	{
		for (int i = 0; i < N_cv; i++)
		{
			Cv[i] = 0;
		}
	}

	//----------------- AGGIORNA LA FUNZIONE DI CORRELAZIONE

	void update_cv() 
	{
		vcm[0] = 0, vcm[1] = 0, vcm[2] = 0;
		for (int i = 0; i < N; i++)
		{
			vcm[0] += v[i][0];
			vcm[1] += v[i][1];
			vcm[2] += v[i][2];
		}
		vcm[0] /= N;
		vcm[1] /= N;
		vcm[2] /= N;
		int index = Tstat% N_cv;
		if (index == 0)
		{
			
			for (int i = 0; i < N; i++)
			{

				v0[i][0] = v[i][0]-vcm[0];
				v0[i][1] = v[i][1]-vcm[1];
				v0[i][2] = v[i][2]-vcm[2];
			}
		}
		for (int i = 0; i < N; i++)
		{
			Cv[index] += v0[i][0] * (v[i][0]-vcm[0]) + v0[i][1] * (v[i][1] - vcm[1]) + v0[i][2] * (v[i][2] - vcm[2]);
		}
		
	}

	//----------------- CREA IL RETICOLO

	void init_r() 
	{
		double dL = L / n;
		int i = 0;
		double dr[3],f[3], dr2, dr_mod, f_mod, der_pot;
		E[1] = 0,viriale = 0;
		double BCC[2][3] = {{0,0,0},{0.5,0.5,0.5} };
		double FCC[4][3] = {{0,0,0},{0.5,0.5,0},{0.5,0,0.5},{0,0.5,0.5}};
		for (int x = 0; x < n; x++)
		{
			for (int y = 0; y < n; y++)
			{
				for (int z = 0; z < n; z++)
				{
					if (M == 1)
					{
						r[i][0] = -L2 + dL * x;
						r[i][1] = -L2 + dL * y;
						r[i][2] = -L2 + dL * z;
						i++;
					}
					if (M == 2)
					{
						for (int j = 0; j < M; j++)
						{
							r[i][0] = -L2 + dL * ((double)x+BCC[j][0]);
							r[i][1] = -L2 + dL * ((double)y+BCC[j][1]);
							r[i][2] = -L2 + dL * ((double)z+BCC[j][2]);
							i++;
							
						}
					}
					if (M == 4)
					{
						for (int j = 0; j < M; j++)
						{
							r[i][0] = -L2 + dL * ((double)x + FCC[j][0]);
							r[i][1] = -L2 + dL * ((double)y + FCC[j][1]);
							r[i][2] = -L2 + dL * ((double)z + FCC[j][2]);
							
							
							i++;
						}
					}
				}
			}
		}
		for (int i = 0; i < N; i++)
		{
			r[i][0] = cpc(r[i][0] + v[i][0] * dt);
			r[i][1] = cpc(r[i][1] + v[i][1] * dt);
			r[i][2] = cpc(r[i][2] + v[i][2] * dt);
			F[i][0] = 0, F[i][1] = 0, F[i][2] = 0;
		}
		for (int i = 1; i < N; i++)
		{
			for (int j = 0; j < i; j++)
			{
				dr[0] = cpc(r[i][0] - r[j][0]);
				dr[1] = cpc(r[i][1] - r[j][1]);
				dr[2] = cpc(r[i][2] - r[j][2]);
				dr2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
				dr_mod = sqrt(dr2);
				if (dr_mod < L2) {
					der_pot = DLJ(dr_mod);
					f_mod = -der_pot / dr_mod;
					f[0] = f_mod * dr[0];
					f[1] = f_mod * dr[1];
					f[2] = f_mod * dr[2];
					F[i][0] += f[0];
					F[j][0] -= f[0];// Forza opposta;
					F[i][1] += f[1];
					F[j][1] -= f[1];
					F[i][2] += f[2];
					F[j][2] -= f[2];
					E[1] += LJ(dr2);
					viriale -= der_pot * dr_mod;
				}
			
			}
		}
		E[2] = E[0] + E[1];
		Ts = 2 * E[0] / (3 * N);
		pres_i = rho * T + viriale / (3 * V);
	}

	//----------------- CREA LA DISTRIBUZIONE DELLE VELOCITÀ

	void init_v() 
	{
		E[0] = 0;
		double sig = sqrt(T);
		double x1, x2, x3, x4;
		for (int i = 0; i < N; i++)
		{
			x1 = (double)rand() /(RAND_MAX+1.); 
			x2 = (double)rand() /(RAND_MAX+1.); 
			x3 = (double)rand() /(RAND_MAX+1.); 
			x4 = (double)rand() /(RAND_MAX+1.); 
			v[i][0] = sig*sqrt(-2 * log(1 - x1)) * cos(2 * pi * x2);
			v[i][1] = sig*sqrt(-2 * log(1 - x2)) * sin(2 * pi * x1);
			v[i][2] = sig*sqrt(-2 * log(1 - x3)) * sin(2 * pi * x4);
			E[0] += v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
		}
		E[0] *= 0.5;
	}

	//----------------- CONDIZIONE PERIODICA AL CONTORNO

	double cpc(double x) 
	{
		return x - L * rint(x / L);
	}

	//----------------- CALCOLA LE FORZE AGENTI SU OGNI PARTICELLA

	void calcF() 
	{
		int idx;
		double dr[3],f[3], dr2, dr_mod, f_mod, der_pot;
		E[1] = 0;
		viriale = 0;
		for (int i = 0; i < N; i++)
		{
			F[i][0] = 0, F[i][1] = 0, F[i][2] = 0;
		}
		for (int i = 1; i < N; i++)
		{
			for (int j = 0; j < i; j++)
			{
				dr[0] = cpc(r[i][0] - r[j][0]);
				dr[1] = cpc(r[i][1] - r[j][1]);
				dr[2] = cpc(r[i][2] - r[j][2]);
				dr2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
				dr_mod = sqrt(dr2);
				if (dr_mod < L2) {
					der_pot = DLJ(dr_mod);
					f_mod = -der_pot / dr_mod;
					f[0] = f_mod * dr[0], f[1] = f_mod * dr[1], f[2] = f_mod * dr[2];
					F[i][0] += f[0];
					F[j][0] -= f[0];// Forza opposta;
					F[i][1] += f[1];
					F[j][1] -= f[1];
					F[i][2] += f[2];
					F[j][2] -= f[2];
					E[1] += LJ(dr2);
					viriale -= der_pot * dr_mod;
					if (fase_stat)
					{
						idx = (int)(dr_mod / dR);
						gr[idx] += 2;
					}
				}
			}
		}
	}

	//----------------- AGGIORNA LE POSIZIONI DELLE PARTICELLE 
	//					TRAMITE L'ALGORITMO VELOCITY-VERLET

	void passo_verlet() 
	{
			E[0] = 0;
			double dt2 = dt * 0.5;
			double dt22 = dt2 * dt;
			for (int i = 0; i < N; i++)
			{
				r[i][0] = cpc(r[i][0] + v[i][0] * dt + F[i][0] * dt22);
				r[i][1] = cpc(r[i][1] + v[i][1] * dt + F[i][1] * dt22);
				r[i][2] = cpc(r[i][2] + v[i][2] * dt + F[i][2] * dt22);
				v[i][0] += F[i][0] * dt2;
				v[i][1] += F[i][1] * dt2;
				v[i][2] += F[i][2] * dt2;
			}
			calcF();
			for (int i = 0; i < N; i++)
			{
				v[i][0] += F[i][0] * dt2;
				v[i][1] += F[i][1] * dt2;
				v[i][2] += F[i][2] * dt2;
				E[0] += v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
			}
			E[0] *= 0.5;
			E[2] = E[0] + E[1];
			Ts = 2 * E[0] / (3 * N);
			pres_i = rho * T + viriale / (3 *V);
		    t += dt;
	}

	//----------------- CREA IL VETTORE DI DISTRIBUZIONE DI COPPIA

	void init_Gr() 
	{
		for (int i = 0; i < N_bin; i++)
		{
			norm[i] = 4. / 3 * pi * (pow((i + 1) * dR, 3) - pow(i * dR, 3)) * rho * N;
			gr[i] = 0;
		}
	}
};


