import numpy as np
from qutip import *

class Hamiltonian:
    '''
    Class including most important info about the hamiltonian
    
    Args:

        wa (float) : atom frequency
        wc (float) : cavity frequency
        Omega (float) : Rabi oscillation frequency
        w (float) : driving frequency on atom
        f (float) : amplitude of driving on atom
        N (float) : numver of photon states (# of photons + 1)
        rwa (bool) : Opt (default = False), Employ RWA?
        driven (bool) : Opt (default = False), drive atom?
        H_JC (QObj) : Hamiltonian generated

    Attributes:
        _Hjc (N, rwa, driven)
    '''
    def __init__(self, wa, wc, Omega, w, f, N, rwa = False, driven = False, rwa_drive = False):
        self.wa = wa
        self.wc = wc
        self.Omega = Omega
        self.w = w
        self.f = f
        self._Hjc(N, rwa = rwa, driven =driven, rwa_drive = rwa_drive)
    def _Hjc (self, N, rwa = False, driven = False, rwa_drive = False):
        # N = number of photon states (number of photons + 1(no photon state) )
        Omega = self.Omega
        wc = self.wc
        wa = self.wa
        w = self.w
        f = self.f
        sigma_z = - tensor(qeye(N), sigmaz())
        a = tensor(destroy(N), qeye(2))
        a_dag = a.dag()
        sm = tensor(qeye(N),destroy(2))
        sigma_x = tensor(qeye(N),sigmax()) #for the qubit system
        sigma_p = tensor(qeye(N),sigmam())
        sigma_m = tensor(qeye(N),sigmap())

        #either that or
        #psi0 = tensor(fock(N, 0), fock(2,1) )
        #psi1 = tensor(fock(N, 1), fock(2,1) )   
        H0 = wc * a.dag() * a + wa/2 * sigma_z
        coeff = coefficient(lambda t: np.cos(w*t))
        coeff_1 = coefficient(lambda t: np.exp(-1j*w*t))
        coeff_2 = coefficient(lambda t: np.exp(1j*w*t))


        if rwa == True:
            if driven == True:
                if rwa_drive == False:
                    self.H_JC = (Omega/2*((a.dag()*sm) + (sm.dag()*a)) + H0 + (f*sigma_x)*coeff)
                else:
                    self.H_JC = (Omega/2*((a.dag()*sm) + (sm.dag()*a)) + H0 + (f/2)*(sigma_p*coeff_1 + sigma_m*coeff_2))
            else:
                self.H_JC = (Omega/2*((a.dag()*sm) + (sm.dag()*a)) + H0)
        else:
            if driven == True:
                if rwa_drive == False:
                    self.H_JC = (Omega/2*(a.dag()+a) * (sm+sm.dag()) + H0 + (f*sigma_x)*coeff)
                else:
                    self.H_JC = (Omega/2*(a.dag()+a) * (sm+sm.dag()) + H0 + (f/2)*(sigma_p*coeff_2 + sigma_m*coeff_1))
            else:
                self.H_JC = (Omega/2*(a.dag()+a) * (sm+sm.dag()) + H0)