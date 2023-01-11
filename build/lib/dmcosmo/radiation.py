"""
Contains luminosity function calculation, stellar parameters
"""

from .constants import *
from .massfct import *
from .cosmo import H_hubble

nu_min, nu_max =nu_al,  nu_LL
M_UV_max= -15

def S_fct(Mh, Mt, g3, g4):
    return (1 + (Mt / Mh) ** g3) ** g4


def f_star_Halo(param):
    f_st = param.rad.f_st
    Mp = param.rad.Mp
    g1 = param.rad.g1
    g2 = param.rad.g2
    Mt = param.rad.Mt
    g3 = param.rad.g3
    g4 = param.rad.g4
    return 2 * f_st / ((Mh / Mp) ** g1 + (Mh / Mp) ** g2) * S_fct(Mh, Mt, g3, g4)


class Radiation():
    """
    M_st_dot  : MAR or tstar (dMst/dt = fst*dMh/dt or dMst/dt = Mst*H(t)/tstar)
    """

    def __init__(self, param):
        self.z = param.code.z
        self.N_al = param.rad.N_al

        HMF = HMF(param)
        HMF.generate_HMF(param)

        self.tab_Mh = HMF.tab_M
        self.MassFct = HMF.HMF

        self.f_star = f_star_Halo(param)

        if not (self.f_star <= 1).all():
            print('warning, some fstar values are larger than one.')

        self.M_st_dot = param.code.Maccr
        self.tstar = 1 # param.code.tstar
        self.K_UV = self.K__UV()
        self.M_UV_max = M_UV_max
        self.alph_mdot = param.code.alph_mdot
        self.H = H_hubble(self.z,param)


    def fstar(self):  ### careful, Ob/Om not included in fstar
        return self.f_star

    def K__UV(self):  ### (Msol/yr) / (erg/s/Hz)    ### WE HAVE TO CHANGE THAT TO THE PROPER PARAMETRIZED STELLAR/QASAR SPECTRUM
        alS = 0.86    ### power law lyman alpha SED
        Anorm = (1 - alS) / (nu_LL ** (1 - alS) - nu_al ** (1 - alS))
        Enu = 1 / (nu_max - nu_min) * integrate.quad(lambda nu: hP * Anorm * nu ** (-alS + 1), nu_min, nu_max)[0]  # mean UV energy
        return m_pr * sec_per_yr / (kg_per_Msol * self.N_al * (Enu / eV_per_erg))

    def Luminosity_(self):

        alpha = self.alph_mdot  ##for exp mass accretion
        dM_dz = alpha * self.tab_Mh  ### Exp Accr
        dM_dt = dM_dz * self.H * (self.z + 1)
        dMstar_dz = dM_dz * self.f_star * Ob / Om
        if self.M_st_dot == 'MAR':
            dMstar_dt = dM_dt * self.f_star * Ob / Om

        if self.M_st_dot == 'tstar':
            dMstar_dt = self.tab_Mh * self.f_star * Ob / Om * H_hubble(self.z) / self.tstar

        L_UV = dMstar_dt / self.K_UV  ## L_UV in erg/s/Hz
        M_UV = 51.63 - 0.4 ** -1 * np.log10(L_UV)
        dMstardot_dMh = np.gradient(dMstar_dz, self.tab_Mh)
        dMuv_dLuv = 1 / (0.4 * np.log(10) * dMstar_dz)
        dMUV_dMh = dMstardot_dMh * dMuv_dLuv
        phi_Muv = h0 ** 3 * self.HMF_ / self.tab_Mh * dMUV_dMh ** -1
        return M_UV, phi_Muv, L_UV  ### magnitude and mag^-1 (Mpc)^-3

    def M_star_dot(self, z_array):
        alpha = 0.79  ##for exp mass accretion
        dM_dz = alpha * self.tab_Mh  ### Exp Accr
        dM_dt = dM_dz * H_hubble(z_array[:, None]) * (z_array[:, None] + 1)
        dMstar_dz = dM_dz * self.f_star * Ob / Om
        dMstar_dt_MAR = dM_dt * self.f_star * Ob / Om
        dMstar_dt_tstar = self.tab_Mh * self.f_star * Ob / Om * H_hubble(z_array[:, None]) / self.tstar
        return dMstar_dt_MAR, dMstar_dt_tstar  ### magnitude and mag^-1 (Mpc)^-3

    def SFRD(self):
        alpha = 0.79  ## For exp mass accretion
        dM_dz = alpha * self.tab_Mh  ## Exp Accr
        dM_dt = dM_dz * H_hubble(self.z) * (self.z + 1)  ##
        sfrd_ = np.trapz(self.f_star * Ob / Om * dM_dt * self.HMF_ / self.tab_Mh, self.tab_Mh)
        return sfrd_  ##Msol per year per (Mpc/h)^3

    def Integrated_UV_Lumi(self):  # Cumulated UV light per Hz, compare to fig 11 of 1509.06764
        sfrd_ = self.SFRD()
        return sfrd_ / self.K_UV  #### in erg/sec/Hz/(Mpc/h)^3

    def Integrated_UV_mag(self):  # Here we integrate the Lumi function with a minimum observational limit M_UV_min
        M_UV, phi_Muv, L_UV = self.Luminosity_()
        indices = np.where(M_UV < self.M_UV_max)
        return -np.trapz(phi_Muv[indices] * L_UV[indices], M_UV[indices])  #### in erg/sec/Hz/(Mpc)^3





