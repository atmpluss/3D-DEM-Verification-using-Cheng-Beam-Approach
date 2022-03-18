import numpy as np
import math

def delta_calculation(F, grain_list, bond_dic):
    # global grain_list
    # global bond_dic

    bond = bond_dic[(0,1)]

    L0 = abs(grain_list[-1].pos[0] - grain_list[0].pos[0])
    Ib = 0.25 * math.pi * bond.Rb**4
    xs = np.array([x.pos[0] for x in grain_list])
    # delta_y = F[1]*np.multiply(xs,xs)/(6. * bond.Eb * Ib) * (3. * L0 - xs) + 10*F[1]*xs/(9*bond.Gb * bond.Ab)

    delta_y = np.multiply((F[1]*np.multiply(xs,xs)/(6. * bond.Eb * Ib)),(3. * L0 - xs) )

    delta_b = 0

    return [delta_y, delta_b,L0]