"""
Mass Function Class
"""

import numpy as np
from math import pi
from .cosmo import *

class HMF:

    def __init__(self, param):
        data = np.loadtxt(param.cosmo.ps)
        self.k_lin = data[:, 0]
        self.P_lin = data[:, 1]
        self.tab_M = np.logspace(np.log10(param.code.m_min), np.log10(param.code.m_max), param.code.Mbin, base=10)   # [Msol/h]
        self.tab_R = ((3 * self.tab_M / (4 * param.cosmo.rho_c * param.cosmo.Om * pi)) ** (1. / 3)) / param.PS.c     # [Mpc/h]
        self.z = param.code.z

    def sigma_square(self,param):
        if param.PS.filter == 'tophat':
            self.sigma2 = np.trapz(self.k_lin ** 2 * self.P_lin * W_tophat(self.k_lin * self.tab_R[:,None]) ** 2 / (2 * pi ** 2) , self.k_lin, axis = 1)
            self.dlnsigmdlnR = np.trapz(self.k_lin ** 3 * self.tab_R[:,None] * self.P_lin * W_tophat(self.k_lin * self.tab_R[:, None]) * derivative_W_tophat(self.k_lin * self.tab_R[:, None])/(2 *self.sigma2[:,None] * pi ** 2), self.k_lin, axis=1)
        elif param.PS.filter == 'sharpk':
            self.sigma2 = np.trapz(self.k_lin ** 2 * self.P_lin * W_sharpk(self.k_lin * self.tab_R[:, None]) ** 2 / (2 * pi ** 2), self.k_lin, axis=1)
            self.dlnsigmdlnR = np.interp(1/self.tab_R, self.k_lin, self.P_lin) /4 /np.pi**2 / self.tab_R**3 / self.sigma2
        elif param.PS.filter == 'smoothk':
            self.sigma2 = np.trapz(self.k_lin ** 2 * self.P_lin * W_smooth_k(self.k_lin * self.tab_R[:, None],param) ** 2 / (2 * pi ** 2), self.k_lin, axis=1)
            self.dlnsigmdlnR = np.trapz(self.k_lin ** 3 * self.tab_R[:,None] * self.P_lin * W_smooth_k(self.k_lin * self.tab_R[:, None],param) * derivative_W_smooth(self.k_lin * self.tab_R[:, None],param)/(2 * self.sigma2[:,None] * pi ** 2), self.k_lin, axis=1)
        else :
            print('filter should be tophat, sharpk or smoothk')

    def generate_HMF(self,param):
        self.sigma_square(param)
        self.sigma_z = D(1. / (self.z  + 1), param) * np.sqrt(self.sigma2)
        self.f_ST = crossing_f_ST(self.sigma_z, param)
        self.HMF = self.f_ST * param.cosmo.rho_c * param.cosmo.Om * np.abs(self.dlnsigmdlnR) / 3 / self.tab_M






