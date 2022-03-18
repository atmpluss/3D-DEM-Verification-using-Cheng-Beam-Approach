import numpy as np
import math


def bonded_grains_interaction(i,j, dt, grain_list, bond_dic):
    id1 = grain_list[i].id
    id2 = grain_list[j].id

    posij = grain_list[i].pos - grain_list[j].pos
    cij = np.sqrt(np.dot(posij, posij))
    dn_current = cij - grain_list[i].r - grain_list[j].r
    dn = dn_current - bond_dic[(id1,id2)].dn_init
    nij = posij

    if(np.linalg.norm(nij) != 0):
        nij /= np.linalg.norm(nij)

    contactPoint = cij/2. * nij + grain_list[j].pos
    rci = contactPoint - grain_list[i].pos
    rcj = contactPoint - grain_list[j].pos

    Ui = grain_list[i].v  + np.cross(rci, grain_list[i].an_v)
    Uj = grain_list[j].v  + np.cross(rcj, grain_list[j].an_v)

    relVelij = Ui - Uj
    relVelijN = nij * np.dot(nij, relVelij)
    relVelijT = relVelij - relVelijN

    fnabs = 0
    fn = np.array([0.,0.,0.])
    if (dn >= 0):
        #traction
        fnabs = -1 * bond_dic[(id1, id2)].kn * dn
        fn = nij * fnabs - bond_dic[(id1, id2)].nun * relVelijN
    else:
        #compression
        fnabs = -1 * bond_dic[(id1, id2)].kn * dn
        if(fnabs < 0.):
            fnabs = 0.
        fn = nij * fnabs - bond_dic[(id1, id2)].nun * relVelijN

    Us_pre = bond_dic[(id1, id2)].us
    Us_pre = Us_pre - nij*(np.dot(nij,Us_pre))
    Us = Us_pre + dt * relVelijT
    fTLS = -1 * bond_dic[(id1, id2)].ks * Us -1 * bond_dic[(id1, id2)].nus*relVelijT
    bond_dic[(id1, id2)].us = Us

    a_prime_ij = ((grain_list[i].r - 0.5 * abs(dn)) * (grain_list[j].r - 0.5 * abs(dn))) /((grain_list[i].r - 0.5 * abs(dn))  + (grain_list[j].r - 0.5 * abs(dn)) )
    a_ij = (grain_list[i].r * grain_list[j].r)/ (grain_list[i].r + grain_list[j].r)

    crossLi = np.cross(nij, grain_list[i].an_v)
    crossLj = np.cross(nij, grain_list[j].an_v)

    vijr = -1 * a_prime_ij * (crossLi-crossLj)

    Ur_pre = bond_dic[(id1, id2)].ur - nij * np.dot(nij,bond_dic[(id1, id2)].ur )
    Ur = Ur_pre + dt * vijr

    fr = -1 * bond_dic[(id1, id2)].kr * Ur - bond_dic[(id1, id2)].nur * vijr
    torque = a_ij * np.cross(nij, fr)
    bond_dic[(id1, id2)].ur = Ur

    vijo = a_ij * nij*(np.dot(nij, grain_list[i].an_v) - np.dot(nij, grain_list[j].an_v))

    Uo_pre = bond_dic[(id1, id2)].uo - nij * np.dot(nij, bond_dic[(id1, id2)].uo)
    Uo = Uo_pre + dt * vijo

    fo = -1 * bond_dic[(id1, id2)].ko * Uo - bond_dic[(id1, id2)].nuo * vijo
    torsion = a_ij * fo
    bond_dic[(id1, id2)].ur = Ur


    # rupture = (np.linalg.norm(fTLS)/bond_dic[(id1, id2)].y_s)**2 + (fnabs/bond_dic[(id1, id2)].y_n) + (np.linalg.norm(torque)/bond_dic[(id1, id2)].y_r)**2 + (np.linalg.norm(torsion)/bond_dic[(id1, id2)].y_o)**2 - 1.
    rupture = -1
    if(rupture >= 0 ):
        fn = np.array([0.,0.,0.])
        fTLS = np.array([0.,0.,0.])
        del(bond_dic[(id1,id2)])
        print("the bond is broken!")
        return

    bond_dic[(id1, id2)].fn = fn
    bond_dic[(id1, id2)].fs = fTLS

    grain_list[i].a += fn + fTLS
    grain_list[j].a += (fn + fTLS) * -1.

    crossthi = np.cross(rci, fTLS)
    FFi = -1 * crossthi
    FFi += torque + torsion

    crossthj = np.cross(rcj, -1*fTLS)
    FFj = -1 * crossthj
    FFj += -1 * (torque + torsion)

    grain_list[i].an_a += FFi
    grain_list[j].an_a += FFj



