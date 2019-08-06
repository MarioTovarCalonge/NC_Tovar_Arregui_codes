# This Python file uses the following encoding: utf-8
import numpy as np
import math as m
from scipy.stats import norm
from scipy import special
import scipy.optimize as optim
import matplotlib.pyplot as plt
import os
import sys
import matplotlib.patches as mpatches
import matplotlib
matplotlib.rcParams.update({'errorbar.capsize': 3})

class calculus():

    def S_c(self, x, N_0, beta):
        return N_0*m.exp(-beta*x)

    def L_c(self, x, N_0, beta, r_L, p, q):
        return (beta*(1-p)*N_0)/(r_L + beta*p*q - beta) * ( m.exp(-beta*x) - m.exp(-(r_L+beta*p*q)*x))

    def F_c(self, x, N_0, beta, r, r_L, p, q):
        K = -beta*p*N_0*((r-r_L-beta*p*q-q*beta*(1-p))/((r-beta)*(r-r_L-beta*p*q)))
        y = K*m.exp(-r*x) + (beta*p*N_0*m.exp(-beta*x))/(r-beta) +((beta*p*q*beta*(1-p)*N_0)/(r_L+beta*p*q-beta))*((m.exp(-beta*x)/(r-beta)) - (m.exp(-(r_L+beta*p*q)*x)/(r-r_L-beta*p*q)))
        return y

    def S_v(self, x, N_0, beta, eps_beta, eps_r, eps_p):
        return N_0*m.exp(-(1-eps_beta)*beta*x)

    def L_v(self, x, N_0, beta, r_L, p, q, eps_beta, eps_r, eps_p):
        return ((1-eps_beta)*beta*(1-(1-eps_p)*p)*N_0)/(r_L + (1-eps_beta)*beta*(1-eps_p)*p*q - (1-eps_beta)*beta) * ( m.exp(-(1-eps_beta)*beta*x) - m.exp(-(r_L+(1-eps_beta)*beta*(1-eps_p)*p*q)*x))

    def F_v(self, x, N_0, beta, r, r_L, p, q, eps_beta, eps_r, eps_p):
        K = -(1-eps_beta)*beta*(1-eps_p)*p*N_0*(((1-eps_r)*r-r_L-(1-eps_beta)*beta*(1-eps_p)*p*q-q*(1-eps_beta)*beta*(1-(1-eps_p)*p))/(((1-eps_r)*r-(1-eps_beta)*beta)*((1-eps_r)*r-r_L-(1-eps_beta)*beta*(1-eps_p)*p*q)))
        y = K*m.exp(-(1-eps_r)*r*x) + ((1-eps_beta)*beta*(1-eps_p)*p*N_0*m.exp(-(1-eps_beta)*beta*x))/((1-eps_r)*r-(1-eps_beta)*beta) +(((1-eps_beta)*beta*(1-eps_p)*p*q*(1-eps_beta)*beta*(1-(1-eps_p)*p)*N_0)/(r_L+(1-eps_beta)*beta*(1-eps_p)*p*q-(1-eps_beta)*beta))*((m.exp(-(1-eps_beta)*beta*x)/((1-eps_r)*r-(1-eps_beta)*beta)) - (m.exp(-(r_L+(1-eps_beta)*beta*(1-eps_p)*p*q)*x)/((1-eps_r)*r-r_L-(1-eps_beta)*beta*(1-eps_p)*p*q)))
        return y

    def D_c(self, x0, y0, beta, r_l, r, p, q):
        return y0 - self.S_c(x0, y0, beta) - self.L_c(x0, y0, beta, r_l, p, q) - self.F_c(x0, y0, beta, r, r_l, p, q)

    def D_v(self, x0, y0, beta, r_l, r, p, q, eps_beta, eps_r, eps_p):
        return y0 - self.S_v(x0, y0, beta, eps_beta, eps_r, eps_p) - self.L_v(x0, y0, beta, r_l, p, q, eps_beta, eps_r, eps_p) - self.F_v(x0, y0, beta, r, r_l, p, q, eps_beta, eps_r, eps_p)

    def error_as_usual(self, N_0, t_max, beta, r_l, r, p, q, eff_p_central, rho, rho_up, rho_down, eff_b, eff_b_up, eff_b_down, eff_r, eff_r_up, eff_r_down):
        delta = [0 for i in range(6)]
        delta_pos = 0.0
        delta_neg = 0.0

        data = (rho_up, N_0, t_max, beta, r_l, r, p, q, eff_b, eff_r_up)
        try:
            delta[0] = float(optim.brentq(funS, -100, 10, args = data)) - eff_p_central
        except:
            delta[0] = 0.0

        data = (rho_down, N_0, t_max, beta, r_l, r, p, q, eff_b, eff_r_down)
        try:
            delta[1] = float(optim.brentq(funS, -100, 10, args = data)) - eff_p_central
        except:
            delta[1] = 0.0

        data = (rho, N_0, t_max, beta, r_l, r, p, q, eff_b_up, eff_r)
        try:
            delta[2] = float(optim.brentq(funS, -100, 10, args = data)) - eff_p_central
        except:
            delta[2] = 0.0

        data = (rho, N_0, t_max, beta, r_l, r, p, q, eff_b_down, eff_r)
        try:
            delta[3] = float(optim.brentq(funS, -100, 10, args = data)) - eff_p_central
        except:
            delta[3] = 0.0

        for i in range(len(delta)):
            if (delta[i]>0):
                delta_pos = delta_pos + delta[i]*delta[i]
            else:
                delta_neg = delta_neg + delta[i]*delta[i]

        delta_pos = m.sqrt(delta_pos) + eff_p_central
        delta_neg = eff_p_central - m.sqrt(delta_neg)
        return delta_neg, delta_pos

    def merge_CI(self, vec, vec_low, vec_up, M):
        grid_size = 10000
        pos_grid = [0 for i in range(grid_size)]
        cum_area_grid = [0 for i in range(grid_size)]
        sigma = [0 for i in range(2)]
        #Locate max and min values
        maxi = max(vec_up)
        mini = min(vec_low)
        #grid
        max_grid = m.log(1 - mini)
        try:
            min_grid = m.log(1 - maxi)
        except:
            maxi = 0.99
            min_grid = m.log(1 - maxi)

        delta = (max_grid - min_grid)/grid_size
        for i in range(grid_size):
            pos_grid[i] = min_grid + i*delta
            cum_area_grid[i] = 0
        #Cargamos las áreas en el grid
        for i in range(M):
            try:
                sigma[1] = 0.5*(m.log(1 - vec_low[i]) - m.log(1 - vec[i]))
            except:
                sigma[1] = 1.0
            try:
                sigma[0] = 0.5*(m.log(1 - vec[i]) - m.log(1 - vec_up[i]))
            except:
                sigma[0] = -1.0
            try:
                mu = m.log(1 - vec[i])
            except:
                vec[i] = 0.99
                mu = m.log(1 - vec[i])
            for j in range(grid_size):
                cum_area_grid[j] += 0.5*(1 + special.erf((pos_grid[j] - mu)/(m.sqrt(2)*sigma[pos_grid[j] > mu])))
        for j in range(grid_size):
            cum_area_grid[j] = cum_area_grid[j]/M
        #Localiza el 2.5% a cada lado
        i = 0
        min_index = 0
        while min_index == 0:
            if (cum_area_grid[i] > 0.025):
                min_index = i
            i = i + 1

        max_global = (0.025 - cum_area_grid[min_index-1])*(pos_grid[min_index] - pos_grid[min_index-1])/(cum_area_grid[min_index] - cum_area_grid[min_index-1]) + pos_grid[min_index-1]
        max_global = 1 - m.exp(max_global)

        i = 0
        median_index = 0
        while median_index == 0:
            if (cum_area_grid[i] > 0.5):
                median_index = i
            i = i + 1

        median = (0.5-cum_area_grid[median_index-1])*(pos_grid[median_index]-pos_grid[median_index-1])/(cum_area_grid[median_index]-cum_area_grid[median_index-1])+pos_grid[median_index-1]
        median = 1 - m.exp(median)

        i = grid_size - 1
        max_index = grid_size - 1
        while max_index == (grid_size -1):
            if (cum_area_grid[i] < 0.975):
                max_index = i
            i = i - 1

        min_global = (0.975 - cum_area_grid[max_index])*(pos_grid[max_index + 1] - pos_grid[max_index])/(cum_area_grid[max_index + 1] - cum_area_grid[max_index]) + pos_grid[max_index]
        min_global = 1 - m.exp(min_global)

        return median, min_global, max_global

def funS(eps_p, *ar):
    rho, N_0, T, beta, r_L, r, p, q, eps_beta, eps_r = ar
    #Cálculo de poblaciones de control
    S_c = N_0*m.exp(-beta*T)
    L_c = (beta*(1-p)*N_0)/(r_L + beta*p*q - beta) * ( m.exp(-beta*T) - m.exp(-(r_L+beta*p*q)*T))
    K = -beta*p*N_0*((r-r_L-beta*p*q-q*beta*(1-p))/((r-beta)*(r-r_L-beta*p*q)))
    F_c = K*m.exp(-r*T) + (beta*p*N_0*m.exp(-beta*T))/(r-beta) +((beta*p*q*beta*(1-p)*N_0)/(r_L+beta*p*q-beta))*((m.exp(-beta*T)/(r-beta)) - (m.exp(-(r_L+beta*p*q)*T)/(r-r_L-beta*p*q)))
    #Fin
    #Cálculo de poblaciones de vac
    S_v = N_0*m.exp(-(1-eps_beta)*beta*T)
    L_v = ((1-eps_beta)*beta*(1-(1-eps_p)*p)*N_0)/(r_L + (1-eps_beta)*beta*(1-eps_p)*p*q - (1-eps_beta)*beta) * ( m.exp(-(1-eps_beta)*beta*T) - m.exp(-(r_L+(1-eps_beta)*beta*(1-eps_p)*p*q)*T))
    K = -(1-eps_beta)*beta*(1-eps_p)*p*N_0*(((1-eps_r)*r-r_L-(1-eps_beta)*beta*(1-eps_p)*p*q-q*(1-eps_beta)*beta*(1-(1-eps_p)*p))/(((1-eps_r)*r-(1-eps_beta)*beta)*((1-eps_r)*r-r_L-(1-eps_beta)*beta*(1-eps_p)*p*q)))
    F_v = K*m.exp(-(1-eps_r)*r*T) + ((1-eps_beta)*beta*(1-eps_p)*p*N_0*m.exp(-(1-eps_beta)*beta*T))/((1-eps_r)*r-(1-eps_beta)*beta) +(((1-eps_beta)*beta*(1-eps_p)*p*q*(1-eps_beta)*beta*(1-(1-eps_p)*p)*N_0)/(r_L+(1-eps_beta)*beta*(1-eps_p)*p*q-(1-eps_beta)*beta))*((m.exp(-(1-eps_beta)*beta*T)/((1-eps_r)*r-(1-eps_beta)*beta)) - (m.exp(-(r_L+(1-eps_beta)*beta*(1-eps_p)*p*q)*T)/((1-eps_r)*r-r_L-(1-eps_beta)*beta*(1-eps_p)*p*q)))
    #Fin
    D_c = N_0 - S_c - L_c - F_c
    D_v = N_0 - S_v - L_v - F_v
    return rho*D_c - D_v

def main():
    #Lectura de parámetros del modelo.
    r = 0.972160
    r_l = 0.00075
    beta = float(sys.argv[3])
    p = float(sys.argv[4])
    q = float(sys.argv[5])
    t_max = float(sys.argv[2])
    N_0 = int(sys.argv[1])
    imax = int(sys.argv[6])

    ##########################################

    E_beta = float(sys.argv[7])
    E_p = float(sys.argv[8])
    E_r = float(sys.argv[9])

    ##########################################

    print("\nValidating trial for parameters: \nbeta = {}, p = {}, r = {}, r_l = {}, q = {}, T = {}".format(beta, p, r, r_l, q, t_max));
    print("______________________________________________________________________________");

    #Bloque de arrays necesarios para analizar el trial.

    r_control = [[0 for i in range(imax)] for j in range(3)]
    r_vac = [[0 for i in range(imax)] for j in range(3)]

    eps_b = [0 for i in range(imax)]
    eps_b_sd = [0 for i in range(imax)]
    eps_beta_low = [0 for i in range(imax)]
    eps_beta_up = [0 for i in range(imax)]

    eps_p = [0 for i in range(imax)]
    eps_p_low = [0 for i in range(imax)]
    eps_p_up = [0 for i in range(imax)]
    eps_p_sd = [0 for i in range(imax)]

    eps_r = [0 for i in range(imax)]
    eps_r_low = [0 for i in range(imax)]
    eps_r_up = [0 for i in range(imax)]
    eps_r_sd = [0 for i in range(imax)]

    rho = [[0 for i in range(imax)]for j in range(4)]

    #Fin

    app = calculus()

    #Criterio para deshechar trials sobre eps_r

    def eff_r_calc(t_max, N_0, rho, beta, r_l, r, p, q, eff_beta, eff_p):
        min = 1000
        for eff_r_local in np.arange(0.0, 1.0, 0.01):
            #V = app.Bisection(-1.0, 1.0, 0.0001, t_max, N_0, rho, beta, r_l, r, p, q, eff_beta, eff_r_local)
            data = (rho, N_0, t_max, beta, r_l, r, p, q, eff_beta, eff_r_local)
            #V = float(optim.fsolve(funS, 0.0, args = data))
            try:
                V = float(optim.brentq(funS, -100, 10, args = data))
            except:
                V = 10

            resto = (V - eff_p)*(V - eff_p)

            if(resto < min):
                min = resto
                min_eff_r = eff_r_local

        return min_eff_r


    filename1 = '../IGRA_negative/data/data_{:.2f}_{:.2f}_{:.2f}/eps_reunited/eps_b_reunited.txt'.format(E_beta, E_p, E_r)
    filename2 = '../IGRA_negative/data/data_{:.2f}_{:.2f}_{:.2f}/eps_reunited/eps_r_reunited.txt'.format(E_beta, E_p, E_r)
    filename3 = '../IGRA_negative/data/data_{:.2f}_{:.2f}_{:.2f}/eps_reunited/eps_p_reunited.txt'.format(E_beta, E_p, E_r)

    file1 = open(filename1, 'w')
    file2 = open(filename2, 'w')
    file3 = open(filename3, 'w')

    fp = open("Final_values_{}_{:.2f}_{:.2f}-{:.2f}-{:.2f}.txt".format(N_0, t_max, E_beta, E_p, E_r), 'w')
    fp.write("Parameter central value CI_low_value CI_up_value\n")

    #Lectura de los data del análisis con MLE en las dos cohortes.
    c = 0
    with open("../IGRA_negative/data/data_{:.2f}_{:.2f}_{:.2f}/MLE/mle_control.txt".format(E_beta, E_p, E_r), 'r') as f:
        for line in f:
            line = line.split()
            r_control[0][c] = float(line[0])
            r_control[1][c] = float(line[1])
            r_control[2][c] = float(line[2])
            c = c + 1

    c = 0
    info = [0 for i in range(imax)]
    with open("../IGRA_negative/data/data_{:.2f}_{:.2f}_{:.2f}/MLE/mle_vac.txt".format(E_beta, E_p, E_r), 'r') as f:
        for line in f:
            line = line.split()
            r_vac[0][c] = float(line[0])
            if (round(r_vac[0][c],2) == 0.2):
                info[c] = 1
            r_vac[1][c] = float(line[1])
            r_vac[2][c] = float(line[2])
            c = c + 1
    #Fin

    #Lectura del análisis de supervivencia de COX
    c = 0
    with open("../IGRA_negative/data/data_{:.2f}_{:.2f}_{:.2f}/eps_b/cox_eps_beta.txt".format(E_beta, E_p, E_r), 'r') as f:
        for line in f:
            line = line.split()
            eps_b[c] = float(line[0])
            eps_b_sd[c] = abs((float(line[2]) - float(line[1]))/4)
            eps_beta_low[c] = float(line[1])
            eps_beta_up[c] = float(line[2])
            c = c + 1
    #Fin

    #Lectura de la fracción de infectados de un trial respecto al otro
    c = 0
    with open("../IGRA_negative/data/data_{:.2f}_{:.2f}_{:.2f}/rho/rho.txt".format(E_beta, E_p, E_r), 'r') as f:
        for line in f:
            line = line.split()
            rho[0][c] = float(line[0])
            rho[1][c] = float(line[1])
            rho[2][c] = float(line[2])
            rho[3][c] = float(line[3])/2
            c = c + 1
    #Fin

    min_eff_p = 1.0 - 1.0/(1*p)

    eps_b_g = []
    eps_b_g_up = []
    eps_b_g_low = []

    eps_r_g = []
    eps_r_g_up = []
    eps_r_g_low = []

    eps_p_g = []
    eps_p_g_up = []
    eps_p_g_low = []

    rho_g = []
    rho_g_up = []
    rho_g_low = []

    itergood = 0
    delta = [0 for i in range(2)]
    vac_eff_iters = [0 for i in range(4)]
    for i in range(0, imax):
        bias = 0
        eps_r[i] = 1.0 - r_vac[0][i]/r_control[0][i]
        vac_eff_iters[0] = 1.0 - r_vac[1][i]/r_control[0][i]
        vac_eff_iters[1] = 1.0 - r_vac[2][i]/r_control[0][i]
        vac_eff_iters[2] = 1.0 - r_vac[0][i]/r_control[1][i]
        vac_eff_iters[3] = 1.0 - r_vac[0][i]/r_control[2][i]

        delta[0] = 0
        delta[1] = 0
        for j in range(len(vac_eff_iters)):
            delta[(vac_eff_iters[j] - eps_r[i])>0] += (vac_eff_iters[j] - eps_r[i])*(vac_eff_iters[j] - eps_r[i])

        eps_r_low[i] = eps_r[i] - m.sqrt(delta[0])
        eps_r_up[i] = eps_r[i] + m.sqrt(delta[1])

        max_eff_r = eff_r_calc(t_max, N_0, rho[2][i], beta, r_l, r, p, q, eps_b[i], min_eff_p)
        if (eps_r_up[i]>max_eff_r):
            eps_r_up[i] = max_eff_r


        print('Trial {}/{} done'.format(i, imax))
        valid = True
        if (eps_r[i] > max_eff_r):
            #eps_r[i] = max_eff_r
            bias = 1
            valid = False



        data = (rho[0][i], N_0, t_max, beta, r_l, r, p, q, eps_b[i], eps_r[i])
        try:
            eps_p[i] = float(optim.brentq(funS, -100, 10, args = data))
        except:
            eps_p[i] = 0.0
            valid = False
            bias = 6

        eps_p_low[i], eps_p_up[i] = app.error_as_usual(N_0, t_max, beta, r_l, r, p, q, eps_p[i], rho[0][i], rho[2][i], rho[1][i], eps_b[i], eps_beta_up[i], eps_beta_low[i], eps_r[i], eps_r_up[i], eps_r_low[i])
        if (eps_p_low[i] < min_eff_p):
            eps_p_low[i] = min_eff_p

        if (eps_p[i] < min_eff_p):
            eps_p[i] = min_eff_p
            bias = 2
            valid = False
        if (eps_b[i] > 1.0 or eps_p[i]>1):
            bias = 3
            valid = False

        if (abs(eps_r[i]) > 500):
            bias = 4
            valid = False

        if (abs(eps_b[i]) > 500 or abs(eps_p[i]) > 500):
            bias = 5
            valid = False

        if (bias == 0):
            itergood = itergood + 1
            eps_b_g.append(eps_b[i])
            eps_b_g_up.append(eps_beta_up[i])
            eps_b_g_low.append(eps_beta_low[i])

            eps_r_g.append(eps_r[i])
            eps_r_g_up.append(eps_r_up[i])
            eps_r_g_low.append(eps_r_low[i])

            eps_p_g.append(eps_p[i])
            eps_p_g_up.append(eps_p_up[i])
            eps_p_g_low.append(eps_p_low[i])

            rho_g.append(rho[0][i])
            rho_g_up.append(rho[2][i])
            rho_g_low.append(rho[1][i])


        file1.write("{}\t{}\t{}\t{}\t{}\n".format(eps_b[i], eps_beta_low[i], eps_beta_up[i], valid, bias))
        file2.write("{}\t{}\t{}\t{}\t{}\n".format(eps_r[i], eps_r_low[i], eps_r_up[i], valid, bias))
        file3.write("{}\t{}\t{}\t{}\t{}\n".format(eps_p[i], eps_p_low[i], eps_p_up[i], valid, bias))

    #################################################

    for i in range(itergood):
        if(eps_p_g[i] > 1.0):
            eps_p_g[i] = 0.999999
        if eps_p_g_up[i] > 1.0:
            eps_p_g_up[i] = 1 - 1E-16
        if(eps_r_g[i] > 1.0):
            eps_r_g[i] = 0.999999
        if eps_r_g_up[i] > 1.0:
            eps_r_g_up[i] = 1 - 1E-16
    #-----------------------------------------------#
    lower_error = []
    upper_error = []

    mu, mu_min, mu_max = app.merge_CI(eps_b_g, eps_b_g_low, eps_b_g_up, itergood)
    fp.write('eps_beta {}\t{}\t{}\n'.format(np.median(eps_b), mu_min, mu_max))
    lower_error.append(mu - mu_min)
    upper_error.append(mu_max - mu)

    mu, mu_min, mu_max = app.merge_CI(rho_g, rho_g_low, rho_g_up, itergood)

    mu, mu_min, mu_max = app.merge_CI(eps_r_g, eps_r_g_low, eps_r_g_up, itergood)
    fp.write('eps_r {}\t{}\t{}\n'.format(np.median(eps_r), mu_min, mu_max))
    lower_error.append(mu-mu_min)
    upper_error.append(mu_max-mu)

    mu, mu_min, mu_max = app.merge_CI(eps_p_g, eps_p_g_low, eps_p_g_up, itergood)
    fp.write('eps_p {}\t{}\t{}\n'.format(np.median(eps_p), mu_min, mu_max))
    lower_error.append(mu-mu_min)
    upper_error.append(mu_max-mu)

    #################################################

    file1.close()
    file2.close()
    file3.close()
    fp.close()

    e = [lower_error, upper_error]

    medians = [E_beta, E_r, E_p]
    coordinates = [eps_b_g, eps_r_g, eps_p_g]
    medians_data = [round(np.median(eps_b),4), round(np.median(eps_r),4), round(np.median(eps_p),4)]
    medians_data2 = [round(np.median(eps_b),4), round(np.median(eps_r),4), round(np.median(eps_p),4)]

    def set_axis_style(ax, labels):
        ax.get_xaxis().set_tick_params(direction='out')
        ax.xaxis.set_ticks_position('bottom')
        ax.set_xticks(np.arange(1, len(labels) + 1))
        ax.set_xticklabels(labels, fontsize=14)
        ax.set_xlim(0.25, len(labels) + 0.75)
        ax.set_xlabel('')
        ax.set_ylim([-1.5,1.0])

    fs = 10  # fontsize
    pos = [1, 2, 3]

    fig, axes = plt.subplots(nrows=1, ncols=1)

    violin_parts = axes.violinplot(coordinates, pos, widths=0.5,
                      showmeans=False, showextrema=False, showmedians=False)

    for vp in violin_parts['bodies']:
        vp.set_facecolor('dimgrey')

    axes.scatter(pos, medians, marker='o', color='black', s=30, zorder=2, label='Real value')
    axes.scatter(pos, medians_data, marker="_", color='blue', s=400, label='Inference')
    axes.errorbar(pos, medians_data2, yerr=e, color='blue', ls='none')
    #axes.errorbar(pos[0], medians_data2[0], yerr=prom_b, color='blue', ls='none')
    #axes.errorbar(pos[1], medians_data2[1], yerr=prom_r, color='blue', ls='none')
    #axes.errorbar(pos[2], medians_data2[2], yerr=prom_p, color='blue', ls='none')
    labels = ['$\\varepsilon_{\\beta}$', '$\\varepsilon_r$', '$ \\varepsilon_p$']
    set_axis_style(axes, labels)
    #for i in range(len(pos)):
        #axes.text(pos[i]+0.25, medians_data[i]-0.015, medians_data[i], horizontalalignment='center', size='x-small', color='b', weight='semibold')
        #axes.text(pos[i]-0.25, medians[i]-0.015, medians[i], horizontalalignment='center', size='x-small', color='black', weight='semibold')

    fig.subplots_adjust(hspace=0.4)
    handles, labels = axes.get_legend_handles_labels()

    # manually define a new patch
    patch = mpatches.Patch(color='silver', label='Data')
    # handles is a list, so append manual patch
    handles.append(patch)
    # plot the legend
    plt.legend(handles=handles)
    plt.show()


    if (valid):
        output = 0
    else:
        output = 1

    sys.exit(output)


if __name__ == '__main__':
    main()
