#!/usr/bin/Rscript

args<-commandArgs(TRUE)
print(args[1])
promedio = as.numeric(args[1])
eps_beta = as.numeric(args[2])
eps_p = as.numeric(args[3])
eps_r = as.numeric(args[4])

packages_available=rownames(installed.packages())
if(!("survival" %in% packages_available))
install.packages("survival")
library(survival)

valor <- integer(promedio)
ups <- integer(promedio)
lows <- integer(promedio)



for(i in 1:promedio){
    name_file <- paste0(sprintf("../IGRA_negative/data/data_%.2f_%.2f_%.2f/Times", eps_beta, eps_p, eps_r), sprintf("/discrete_times_%03d.txt", i))

    surv_times=read.table(name_file,header=T)
    attach(surv_times)

    cox <- coxph(Surv(Time_inf, Delta_inf)~Cohort,method="efron")

    cox.val.inf=1-exp(coef(cox))
    cox.sd.inf=cox.val.inf*2*sqrt(cox$var)
    sum<-summary(cox)
    cox.low.inf=1-sum$conf.int[4]
    cox.up.inf=1-sum$conf.int[3]

    ups[i] = cox.up.inf
    lows[i] = cox.low.inf
    valor[i] = cox.val.inf

    detach(surv_times)
}
name_file <- paste0(sprintf("../IGRA_negative/data/data_%.2f_%.2f_%.2f/eps_b", eps_beta, eps_p, eps_r), sprintf("/cox_eps_beta.txt"))
Datos <- data.frame(valor, lows, ups)
write.table(Datos, name_file, row.names=FALSE, col.names=FALSE)
