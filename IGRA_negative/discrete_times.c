/* discrete_times.c */
/* Discretiza los tiempos de inf y dis */
/* (se testea a los individuos con una cierta regularidad)*/


#include "stdio.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define delta_t 0.25 //3 meses

int N, promedio;
double T_max;
double eps_beta, eps_p, eps_r;

int main(int argc, char **argv)
{
	if(argc==7)
	{
		sscanf(argv[1],"%d",&N);
		sscanf(argv[2],"%lf",&T_max);
        sscanf(argv[3],"%d", &promedio);
        eps_beta = atof(argv[4]);
        eps_p = atof(argv[5]);
        eps_r = atof(argv[6]);
	}
	else
	{
		printf("Incorrect number of arguments\n");
		exit(1);
	}

	FILE *fp,*fp2;
	int i,i_aux,delta_inf,delta_dis,cohort,count_t;
	double time_inf,time_dis, time_aux, time2_aux;
	char file_name[1000];


    int m;
    for (m = 1; m < promedio+1; m++){
        sprintf(file_name, "../IGRA_negative/data/data_%.2f_%.2f_%.2f/Times/times_beta_%03d.txt", eps_beta, eps_p, eps_r, m);
        fp=fopen(file_name, "rt");
        sprintf(file_name,"../IGRA_negative/data/data_%.2f_%.2f_%.2f/Times/discrete_times_%03d.txt", eps_beta, eps_p, eps_r, m);
        fp2=fopen(file_name,"wt");

    	fscanf(fp,"Ind Time_inf Delta_inf Time_dis Delta_dis Cohort\n");
    	fprintf(fp2,"Ind Time_inf Delta_inf Time_dis Delta_dis Cohort\n");


    	for(i=0;i<2*N;i++)
    	{
    		fscanf(fp,"%d %lf %d %lf %d %d\n",&i_aux,&time_inf,&delta_inf,&time_dis,&delta_dis,&cohort);

    		if(time_inf>(T_max-0.05))
    			fprintf(fp2,"%d %lf %d %lf %d %d\n", i_aux,T_max,0,T_max,0,cohort);
    		else
    		{
    			if(delta_dis==0)
    			{
    				count_t=(int)(time_inf/delta_t);
    				time_inf=(count_t+0.5)*delta_t;
    				fprintf(fp2,"%d %lf %d %lf %d %d\n", i_aux,time_inf,delta_inf,T_max,delta_dis,cohort);
    			}
    			else
    			{
                    if(time_dis==0.0){
                        count_t=(int)(time_inf/delta_t);
        				time_inf=(count_t+0.5)*delta_t;
        				count_t=(int)(time_dis/delta_t);
        				time_dis=(count_t+0.5)*delta_t;
        				fprintf(fp2,"%d %lf %d %lf %d %d\n", i_aux,time_inf,delta_inf,T_max,0,cohort);
                    }
                    else{
        				count_t=(int)(time_inf/delta_t);
        				time_inf=(count_t+0.5)*delta_t;
        				count_t=(int)(time_dis/delta_t);
        				time_dis=(count_t+0.5)*delta_t;
        				fprintf(fp2,"%d %lf %d %lf %d %d\n", i_aux,time_inf,delta_inf,time_dis,delta_dis,cohort);
                    }
    			}
    		}
    	}
    	fclose(fp);
    	fclose(fp2);


        sprintf(file_name, "../IGRA_negative/data/data_%.2f_%.2f_%.2f/Times/complete_times_%03d-control.txt", eps_beta, eps_p, eps_r, m);
        fp=fopen(file_name, "rt");
        sprintf(file_name,"../IGRA_negative/data/data_%.2f_%.2f_%.2f/Times/discrete_times_control_%03d.txt", eps_beta, eps_p, eps_r, m);
        fp2=fopen(file_name,"wt");

    	fscanf(fp,"Individual Time_inf Time\n");
    	fprintf(fp2,"Individual Time_inf Time\n");


    	for(i=0; feof(fp)==0; i++)
    	{
    		fscanf(fp,"%d %lf %lf\n",&i_aux,&time_inf,&time_dis);
    		time_dis = time_inf+time_dis;

    		count_t=(int)(time_inf/delta_t);
    		time_inf=(count_t+0.5)*delta_t;
    		count_t=(int)(time_dis/delta_t);
    		time_dis=(count_t+0.5)*delta_t;

    		time_dis=time_dis-time_inf;
    		fprintf(fp2,"%d %lf %lf\n",i_aux,time_inf,time_dis);
    	}
    	fclose(fp);
    	fclose(fp2);

        sprintf(file_name, "../IGRA_negative/data/data_%.2f_%.2f_%.2f/Times/complete_times_%03d-vac.txt", eps_beta, eps_p, eps_r, m);
        fp=fopen(file_name, "rt");
        sprintf(file_name,"../IGRA_negative/data/data_%.2f_%.2f_%.2f/Times/discrete_times_vac_%03d.txt", eps_beta, eps_p, eps_r, m);
        fp2=fopen(file_name,"wt");

    	fscanf(fp,"Individual Time_inf Time\n");
    	fprintf(fp2,"Individual Time_inf Time\n");

    	for(i=0;feof(fp)==0;i++)
    	{
    		fscanf(fp,"%d %lf %lf\n",&i_aux,&time_inf,&time_dis);
    		time_dis=time_inf+time_dis;

    		count_t=(int)(time_inf/delta_t);
    		time_inf=(count_t+0.5)*delta_t;
    		count_t=(int)(time_dis/delta_t);
    		time_dis=(count_t+0.5)*delta_t;

    		time_dis=time_dis-time_inf;
    		fprintf(fp2,"%d %lf %lf\n",i_aux,time_inf,time_dis);
    	}
    	fclose(fp);
    	fclose(fp2);
    }
}
