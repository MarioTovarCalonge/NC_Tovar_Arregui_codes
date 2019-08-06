#!/usr/bin/Rscript

args<-commandArgs(TRUE)
print(args[1])
promedio = as.numeric(args[2])
t_max = as.numeric(args[1])
eps_beta = as.numeric(args[3])
eps_p = as.numeric(args[4])
eps_r = as.numeric(args[5])

packages_available=rownames(installed.packages())
if(!("stats4" %in% packages_available))
install.packages("stats4")
if(!("bbmle" %in% packages_available))
install.packages("bbmle")

#mle.R

library(bbmle)
library(stats4)

T=t_max

valor <- integer(promedio)
ups <- integer(promedio)
lows <- integer(promedio)

valor_vac <- integer(promedio)
ups_vac <- integer(promedio)
lows_vac <- integer(promedio)

Poisson<-function(x,lambda,T_limit){
    x_prime=x+(1/lambda)-(T_limit*exp(-lambda*T_limit))/(1-exp(-lambda*T_limit))
   	P=lambda*exp(-lambda*x_prime)/(1-exp(-lambda*T_limit))
   	return(P)
}

LL<-function(lambda){
    R=Time-(1/lambda)+((T-Time_inf)*exp(-lambda*(T-Time_inf)))/(1-exp(-lambda*(T-Time_inf)))
    R=Poisson(R,lambda,(T-Time_inf))
    -sum(log(R))
}

for(i in 1:promedio){
    name_file <- paste0(sprintf("../IGRA_negative/data/data_%.2f_%.2f_%.2f/Times", eps_beta, eps_p, eps_r), sprintf("/discrete_times_control_%03d.txt", i))

    surv_times=read.table(name_file, header=T)
    attach(surv_times)

    fit<-mle2(LL,start=list(lambda=1),method="L-BFGS-B",lower=c(lambda=0.000001))
     
    lambda=coef(fit)
    lambda.ci=confint(fit)
    lambda.low=lambda.ci[1]
    lambda.up=lambda.ci[2]
    bias=0

    if(is.na(lambda.low)) {lambda.low=0.01}


    ups[i] = lambda.up
    lows[i] = lambda.low
    valor[i] = lambda

    detach(surv_times)

    name_file <- paste0(sprintf("../IGRA_negative/data/data_%.2f_%.2f_%.2f/Times", eps_beta, eps_p, eps_r), sprintf("/discrete_times_vac_%03d.txt", i))

    surv_times=read.table(name_file, header=T)
    attach(surv_times)
    fit<-mle2(LL,start=list(lambda=0),method="L-BFGS-B",lower=c(lambda=0.02))

    lambda=coef(fit)
    lambda.ci=confint(fit)
    lambda.low=lambda.ci[1]
    lambda.up=lambda.ci[2]

    if(is.na(lambda.low)) {lambda.low=0.01}


    ups_vac[i] = lambda.up
    lows_vac[i] = lambda.low
    valor_vac[i] = lambda

    detach(surv_times)

}

Datos_control <- data.frame(valor, lows, ups)
name_file <- paste0(sprintf("../IGRA_negative/data/data_%.2f_%.2f_%.2f/MLE", eps_beta, eps_p, eps_r), sprintf("/mle_control.txt"))
write.table(Datos_control, name_file, row.names=FALSE,col.names=FALSE)

Datos_vac <- data.frame(valor_vac, lows_vac, ups_vac)
name_file <- paste0(sprintf("../IGRA_negative/data/data_%.2f_%.2f_%.2f/MLE", eps_beta, eps_p, eps_r), sprintf("/mle_vac.txt"))
write.table(Datos_vac, name_file, row.names=FALSE,col.names=FALSE)




