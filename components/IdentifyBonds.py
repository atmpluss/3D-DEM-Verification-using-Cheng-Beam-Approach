import numpy as np
import math
from .BondClass import Bond


def identify_cohesive_bonds(grain_array, E, R, rhoS, en, bond_dic):
    'this function identifies the cohesion bonds at the first time step'

    # identifying cohesive bonds:
    # global bond_dic

    for i in range(len(grain_array)):
        for j in range(i + 1, len(grain_array)):

            posij = grain_array[i].pos - grain_array[j].pos
            cij = np.sqrt(np.dot(posij,posij))
            dn = cij - grain_array[i].r - grain_array[j].r


            if (dn <= 0.01):

                m_eff = grain_array[i].mass() * grain_array[j].mass() / (grain_array[i].mass() + grain_array[j].mass())

                poisson = 0.2
                Eb = E
                Gb = Eb/(2*(1+poisson))
                Rb = R
                Ab = math.pi * Rb**2
                Lb =  cij
                Ib = 0.25 * math.pi * Rb**4
                phi = 20. * (Rb**2) * (1. + poisson)/(3. * Lb**2)

                print(f"hardcoging lb: {0.4} and rb: {Rb}")
                print("id1: " + str(i), " id2: " + str(j) + " dn: " + str(dn))
                Lb = 0.4
                kappa = 2.*(1. + poisson)
                kn = Eb * Ab/Lb
                # ks = kn/kappa
                ks = 12 * Eb * Ib /(Lb**3)

                kr = kn
                ko = kn

                damp = -1 * np.log(en) / np.sqrt(np.log(en) * np.log(en) + math.pi * math.pi)
                dampN = 2 * np.sqrt(kn * m_eff) * damp

                nun = dampN
                nus = nun * 1
                nur = nun * 0.5
                nuo = nun * 0.5

                b0 = Bond(gi=grain_array[i],gj=grain_array[j],initial_dn=dn,dn=0,kn=kn,ks=ks,kr=kr,ko=ko,nu=nun,nus=nus,nur=nur,nuo=nuo, Rb=Rb, Ab=Ab, Ib=Ib,Eb=Eb,Gb=Gb,phi=phi)
                bond_dic[(grain_array[i].id,grain_array[j].id)] = b0



