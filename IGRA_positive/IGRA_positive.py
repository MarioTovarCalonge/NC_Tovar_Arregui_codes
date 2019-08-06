# This Python file uses the following encoding: utf-8
import numpy as np
import math as m
from scipy.stats import norm
from scipy import special
import scipy.optimize as optim
import matplotlib.pyplot as plt
import sys
import matplotlib.patches as mpatches
import matplotlib
matplotlib.rcParams.update({'errorbar.capsize': 3})


class calculus():
    def L_c(self, x, N_0, beta, r_L, p, q):
        return (beta*(1-p)*N_0)/(r_L + beta*p*q - beta) * ( m.exp(-beta*x) - m.exp(-(r_L+beta*p*q)*x))

    def F_c(self, x, N_0, beta, r, r_L, p, q):
        K = -beta*p*N_0*((r-r_L-beta*p*q-q*beta*(1-p))/((r-beta)*(r-r_L-beta*p*q)))
        y = K*m.exp(-r*x) + (beta*p*N_0*m.exp(-beta*x))/(r-beta) +((beta*p*q*beta*(1-p)*N_0)/(r_L+beta*p*q-beta))*((m.exp(-beta*x)/(r-beta)) - (m.exp(-(r_L+beta*p*q)*x)/(r-r_L-beta*p*q)))
        return y

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

    def Bisection(self, a, b, tol, T, y0, rho, beta, r_l, r, p, q, eps_beta, eps_r):
    	c = (a+b)/2.0
    	while (b-a)/2.0 > tol:
    		if self.f(rho, y0, T, beta, r_l, r, p, q, eps_beta, eps_r, c) == 0:
    			return c
    		elif self.f(rho, y0, T, beta, r_l, r, p, q, eps_beta, eps_r, a)*self.f(rho, y0, T, beta, r_l, r, p, q, eps_beta, eps_r, c) < 0:
    			b = c
    		else :
    			a = c
    		c = (a+b)/2.0
    	return c

def funS(eps_p, *ar):
    rho, N_0, L_0, F_0, T, beta, r_L, r, p, q, eps_r = ar

    a_c = (beta*p*q + r_L)
    a_v = (beta*p*q*(1-eps_p) + r_L)
    den_c = r - beta*p*q - r_L
    den_v = (1-eps_r)*r - beta*p*q*(1-eps_p) - r_L

    L_c = L_0*m.exp(-a_c*T)
    L_v = L_0*m.exp(-a_v*T)
    F_c = F_0*m.exp(-r*T) - beta*p*q*L_0*m.exp(-r*T)/den_c + beta*p*q*L_0*m.exp(-a_c*T)/den_c
    F_v = F_0*m.exp(-r*(1-eps_r)*T) - beta*p*q*(1-eps_p)*L_0*m.exp(-r*(1-eps_r)*T)/den_v + beta*p*q*(1-eps_p)*L_0*m.exp(-a_v*T)/den_v

    D_c = N_0 - L_c - F_c
    D_v = N_0 - L_v - F_v

    #D_c = N_0 - L_0*m.exp(-a_c*T) - F_0*m.exp(-r*T) - beta*p*q*L_0*m.exp(-r*T)/(beta*p*q+r_L-r) + beta*p*q*L_0*m.exp(-a_c*T)/(beta*p*q + r_L-r)
    #D_v = N_0 - L_0*m.exp(-a_v*T) - F_0*m.exp(-r*(1-eps_r)*T) - beta*p*q*(1-eps_p)*L_0*m.exp(-r*(1-eps_r)*T)/den_v + beta*p*q*(1-eps_p)*L_0*m.exp(-a_v*T)/den_v
    return rho*D_c - D_v


def main():
    app = calculus()
    t_max = 4
    N_0 = 3000
    eps_b = 0.00
    beta = 0.069
    p = 0.15
    q = 0.21 #Pg180 0.21 Comprobar old=0.65
    r_l = 0.00075
    r = 0.97
    VE_dis = 0.50
    rho = 1.0-VE_dis


    for c in range(8):
        if (c==0):
            mult = 1.0 #0.987
        if (c==1):
            mult = 0.995 #0.987
        if (c==2):
            mult = 0.99 #0.987
        if (c==3):
            mult = 0.987 #0.987
        if (c==4):
            mult = 0.98 #0.987
        if (c==5):
            mult = 0.95 #0.987
        if (c==6):
            mult = 0.80 #0.987
        if (c==7):
            mult = 0.01 #0.987

        L_0 = N_0*(mult)
        F_0 = N_0 - L_0

        filename = 'curve_L{:.2f}_F{:.2f}.txt'.format(100*(mult), 100*(1-mult))
        x = []
        y = []

        f = open(filename, 'w')

        for i in np.arange(0.0, 1.0, 0.001):
            x.append(i)
            data = (rho, N_0, L_0, F_0, t_max, beta, r_l, r, p, q, i)
            eps_p = float(optim.fsolve(funS, 0.0, args = data))
            y.append(eps_p)
            f.write('{} {}\n'.format(i, eps_p))


        if(c==7):
            plt.plot(x, y, label='Balance L={:.1f}% F={:.1f}%'.format(0, 100))
        else:
            plt.plot(x, y, label='Balance L={:.1f}% F={:.1f}%'.format(100*(mult), 100*(1-mult)))
        plt.axis([-0.01, 1, -0.01, 1])
        plt.legend()
        plt.grid(True, lw = 1, ls = '--', c = '.75')
        plt.xlabel('$\\varepsilon_{r}$', fontsize=14)
        plt.ylabel('$\\varepsilon_{p}$', fontsize=14)
        plt.xticks(fontsize=11)
        plt.yticks(fontsize=11)


    plt.savefig('Igra_positive_VEdis.png')
    plt.savefig("Igra_positive_VEdis.pdf")
    #plt.show()



    return 0


if __name__ == '__main__':
    main()
