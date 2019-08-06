#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>

#define NormRANu (2.3283063671E-10F)
#define Z_alpha 1.96
//Para pasar una estructura de este tipo en una función hacemos struct Edge *list[] si es declarada como struct Edge *list[].
//En caso de tener un único vector, definimos la estructura como struct Edge *list y se pasa como struct edge *list.

void chemin(int Ntot, int stat[Ntot][4], int *f1, int *f2, int *f3);
int Trial(int iter, int Ntot, float T, double beta, double r_l, double r, double p, double q, bool vaccine_trial, int *flux, double eff1, double eff2, double eff3);

double generador(void);
void iniciar(void);

unsigned int rueda[256];
unsigned char ind_ran,ig1,ig2,ig3;
unsigned int alea;

int main(int argc, char *argv[]){
    int i;

    double r = 0.972160;
    double r_L = 0.00075;
    double beta = atof(argv[3]);
    double p = atof(argv[4]);
    double q = atof(argv[5]);

    float t_max = atof(argv[2]);
    int N_max = atoi(argv[1]);

    char name[2000];

    bool vaccine = true;

    double eps_beta, eps_r, eps_p, beta_vac, p_vac, r_vac;
    int j;
    int imax = atoi(argv[6]);


    eps_beta = atof(argv[7]);
    eps_p = atof(argv[8]);
    eps_r = atof(argv[9]);


    beta_vac = beta*(1.0 - eps_beta);
    p_vac = p*(1.0 - eps_p);
    r_vac = r*(1.0 - eps_r);


    srand(time(NULL));
	iniciar();
    for(i=0; i<1000; i++)
		generador();


    sprintf(name, "../IGRA_negative/data/data_%.02lf_%.02lf_%.02lf/rho/rho.txt", eps_beta, eps_p, eps_r);

    FILE *f;
    f = fopen(name, "wt");
    int D_c = 0;
    int D_v = 0;
    double VE_dis[3];
    double prop_c, prop_v, se;
    int fl[3];
    fl[0] = fl[1] = fl[2] = 0;
    int fl_v[3];
    fl_v[0] = fl_v[1] = fl_v[2] = 0;



    for(i = 1; i<imax+1; i++){

        if (i==1){
            printf("Trial for parameters: beta = %lf, p = %lf, r = %lf, r_l = %lf, q = %lf, T = %.02f\n", beta, p, r, r_L, q, t_max);
            printf("______________________________________________________________________________\n");
        }

        D_c = Trial(i, N_max, t_max, beta, r_L, r, p, q, false, fl, eps_beta, eps_p, eps_r);

        if(vaccine){

            D_v = Trial(i, N_max, t_max, beta_vac, r_L, r_vac, p_vac, q, true, fl_v, eps_beta, eps_p, eps_r);

        }

        VE_dis[0] = D_v/((double) D_c);
    	prop_c = (double)D_c/N_max;
    	prop_v = (double)D_v/N_max;
    	se = sqrt((double)(1-prop_v)/(D_v)+(double)(1-prop_c)/(D_c));
    	VE_dis[2] = (double)D_v/D_c*exp(Z_alpha*se);
    	VE_dis[1] = (double)D_v/D_c*exp(-Z_alpha*se);

        fprintf(f, "%lf \t %lf \t %lf \t %lf\n", VE_dis[0], VE_dis[1], VE_dis[2], se);
        printf("Iteration %03d/%d completed.\n", i, imax);

    }
    fclose(f);
    printf("\nDone!\n");

}

void chemin(int Ntot, int stat[Ntot][4], int *f1, int *f2, int *f3){
    int FluxS_F_D = 0;
    int FluxS_L_F_D = 0;
    int FluxS_L_D = 0;
    int i;

    for (i = 0; i < Ntot; i++){
        if (stat[i][3] == 1){
            if (stat[i][2] == 1 && stat[i][1] == 0)
                FluxS_L_D+=1;
            if (stat[i][2] == 1 && stat[i][1] == 1)
                FluxS_L_F_D+=1;
            if (stat[i][1] == 1 && stat[i][2] == 0)
                FluxS_F_D+=1;
        }
    }
    *f1 = FluxS_F_D;
    *f2 = FluxS_L_F_D;
    *f3 = FluxS_L_D;
}

int Trial(int iter, int Ntot, float T, double beta, double r_l, double r, double p, double q, bool vaccine_trial, int *flux, double eff1, double eff2, double eff3){
    int i, j, S, L, F, D, flux1, flux2, flux3, D_max;
    int node[Ntot];
    int states[Ntot][4];
    for(j=0; j<Ntot; j++){
        node[j]=0;
        states[j][0] = 0;
        states[j][1] = 0;
        states[j][2] = 0;
        states[j][3] = 0;
    }

    int times[Ntot][2];
    int times_beta[Ntot];
    for (i = 0; i < Ntot; i++){
        times[i][0] = 0;
        times[i][1] = 0;
        times_beta[i] = 0;
    }

    int paso = 365;
    int cohort = 0;
    int aux = 0;

    if (vaccine_trial){
        cohort = 1;
        aux = Ntot;
    }

    beta = beta/((double)paso);
    r_l = r_l/((double)paso);
    r = r/((double)paso);

    char filename[10000];
    FILE *z;
    FILE *d;
    if(vaccine_trial){
        sprintf(filename, "../IGRA_negative/data/data_%.02lf_%.02lf_%.02lf/Times/complete_times_%03d-vac.txt", eff1, eff2, eff3,  iter);
        z=fopen(filename, "wt");
        fprintf(z,"Individual Time_inf Time\n");
        sprintf(filename, "../IGRA_negative/data/data_%.02lf_%.02lf_%.02lf/Times/times_beta_%03d.txt", eff1, eff2, eff3, iter);
        d=fopen(filename, "a");
    }
    else{
        sprintf(filename, "../IGRA_negative/data/data_%.02lf_%.02lf_%.02lf/Times/complete_times_%03d-control.txt", eff1, eff2, eff3, iter);
        z=fopen(filename, "wt");
        fprintf(z,"Individual Time_inf Time\n");
        sprintf(filename, "../IGRA_negative/data/data_%.02lf_%.02lf_%.02lf/Times/times_beta_%03d.txt", eff1, eff2, eff3, iter);
        d=fopen(filename, "wt");
        fprintf(d,"Ind Time_inf Delta_inf Time_dis Delta_dis Cohort\n");
    }

    int survey = 90;
    int add = (survey/2);
    add = 0;
    int c = 0;


    for (i = 0; i < paso*T; i++){
        S = L = F = D = 0;
        for (j = 0; j < Ntot; j++){
            if (node[j] == 0){
                if (generador() < beta){
                    if (generador() < p){
                        node[j] = 1;
                    }
                    else{
                        node[j] = 2;
                    }
                }
            }
            if (node[j] == 1){
                if (generador() < r){
                    node[j] = 3;
                }
            }
            if (node[j] == 2){
                if (generador() < beta*q*p){
                    node[j] = 1;
                }
                if (generador() < r_l){
                    node[j] = 3;
                }
            }
            if(states[j][0] == 0 && node[j]==0)
                states[j][0] = 1;
            if(states[j][1] == 0 && node[j]==1)
                states[j][1] = 1;
            if(states[j][2] == 0 && node[j]==2)
                states[j][2] = 1;
            if(states[j][3] == 0 && node[j]==3)
                states[j][3] = 1;

            if (i==c){
                if (node[j]==1 && times[j][0]==0)
                    times[j][0] = i-add;
                if (node[j]==2 && times[j][0]==0)
                    times[j][0] = i-add;
                if (node[j]==1 && times_beta[j]==0)
                    times_beta[j] = i-add;
                if (node[j]==2 && times_beta[j]==0)
                    times_beta[j] = i-add;
                }
            if (i==c){
                if (node[j]==3 && times[j][1]==0)
                    times[j][1] = i-add;
            }

        }
        if(i==c)
            c++;

        for(j = 0; j < Ntot; j++){
            if(node[j] == 0)
                S += 1;
            if(node[j] == 1)
                F += 1;
            if(node[j] == 2)
                L += 1;
            if(node[j] == 3)
                D += 1;
        }

        if (S+F+L+D!=Ntot){
            printf("Error interno, el numero de nodos no se ha conservado.\n");
            break;
        }
    }
    flux1 = flux2 = flux3 = 0;
    chemin(Ntot, states, &flux1, &flux2, &flux3);
    flux[0] = flux1;
    flux[1] = flux2;
    flux[2] = flux3;

    int infect = 0;
    //FILE *a;
    int status = 1;
    for (i = 0; i < Ntot; i++){
        if (times[i][1]!=0 && times[i][0]!=0){
            fprintf(z, "%04d \t %lf \t %lf\n", i, times[i][0]/((double)paso), times[i][1]/((double)paso)-times[i][0]/((double)paso));
        }
        status = 1;
        if (times_beta[i]==0){
            status = 0;
            times_beta[i] = T*paso;
        }
        else{
            infect++;
        }
        fprintf(d,"%04d %lf %d %lf %d %d\n",i + aux, times_beta[i]/((double)paso), (times_beta[i]/((double)paso) < T), times[i][1]/((double)paso), (times[i][1]/((double)paso) < T), cohort);

    }
    fclose(z);
    fclose(d);
    //fclose(a);
    return D;
}

void iniciar(void){
    int i;

    for(i=0; i<256; i++){
		rueda[i]=((rand()<<16)+rand());
	}
    ind_ran=0;
    ig1=0;
    ig2=0;
    ig3=0;
}

double generador(void){
	double r;

	ig1=ind_ran-24;
	ig2=ind_ran-55;
	ig3=ind_ran-61;
	rueda[ind_ran]=rueda[ig1]+rueda[ig2];
	alea=(rueda[ind_ran]^rueda[ig3]);
	ind_ran++;

	r=alea*NormRANu;
	return r;
}
