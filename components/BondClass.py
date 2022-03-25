import numpy as np
import math


class Bond:

    def __init__(self, gi,gj,initial_dn,dn,kn,ks,kr,ko,nu,nus,nur,nuo, Rb, Ab, Ib,Eb,Gb,phi):
        self.graini_id = gi
        self.grainj_id = gj
        self.dn_init = initial_dn
        self.dn = dn
        self.us= np.array([0,0,0])
        self.ur= np.array([0,0,0])
        self.uo = np.array([0,0,0])
        self.f_n = np.array([0,0,0])
        self.f_s = np.array([0,0,0])
        self.torque = np.array([0,0,0])
        self.torsion = np.array([0,0,0])
        self.rupture = -1
        # stiffness
        self.kn = kn
        self.ks = ks
        self.kr = kr
        self.ko = ko
        # damping
        self.nun = nu
        self.nus = nus
        self.nur = nur
        self.nuo = nuo

        # thresholds
        self.y_n = math.inf
        self.y_s = math.inf
        self.y_r = math.inf
        self.y_o = math.inf
        #Chen article
        self.Rb = Rb
        self.Ab = Ab
        self.Ib = Ib
        self.Eb = Eb
        self.Gb = Gb
        self.phi = phi
