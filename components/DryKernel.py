import numpy as np
import math

def dryParticle(i,j,dt,E,mu,en,ks_to_kn,nus_to_nun, grain_list):
    id1 = grain_list[i].id
    id2 = grain_list[j].id

    posij = grain_list[i].pos - grain_list[j].pos
    cij = np.sqrt(np.dot(posij, posij))
    dn = cij - grain_list[i].r - grain_list[j].r
    if dn >= 0:
        return
    kn = 2. * E * (grain_list[i].r * grain_list[j].r) / (grain_list[i].r + grain_list[j].r)
    m_eff = grain_list[i].mass() * grain_list[j].mass() / (grain_list[i].mass() + grain_list[j].mass())

    damp = -1 * np.log(en) / np.sqrt(np.log(en) * np.log(en) + math.pi * math.pi)
    dampN = 2 * np.sqrt(kn * m_eff) * damp
    ks = ks_to_kn  * kn
    dampS = nus_to_nun * dampN
    nij = posij
    if(np.linalg.norm(nij) != 0.):
        nij /= np.linalg.norm(nij)

    contactPoint = cij / 2. * nij + grain_list[j].pos
    rci = contactPoint - grain_list[i].pos
    rcj = contactPoint - grain_list[j].pos

    Ui = grain_list[i].v + np.cross(rci, grain_list[i].an_v)
    Uj = grain_list[j].v + np.cross(rcj, grain_list[j].an_v)

    relVelij = Ui - Uj
    relVelijN = nij * np.dot(nij, relVelij)
    relVelijT = relVelij - relVelijN

    fn = np.array([0., 0., 0.])
    fn = nij * kn * abs(dn) - relVelijN * dampN

    newtangentialDisplacement = dt * relVelijT
    fTLS = -1 * ks * newtangentialDisplacement - dampS * relVelijT

    t = fTLS
    fFrictionAbs = 0.
    if(np.linalg.norm(t) != 0.):
        t /= np.linalg.norm(t)

    fFrictionAbsDynamic = mu * np.linalg.norm(fn)
    if(np.linalg.norm(fTLS) < fFrictionAbsDynamic ):
        fFrictionAbs = fFrictionAbsDynamic
    else:
        fFrictionAbs = fFrictionAbsDynamic
        newtangentialDisplacement = (fFrictionAbs * t - dampS * relVelijT) / ks

    fTabs = min(np.linalg.norm(fTLS), fFrictionAbs)
    fT = fTabs*t

    grain_list[i].a += fn + fT
    grain_list[j].a += -1. * (fn + fT)

    crossthi = np.cross(rci, fT)
    FFi = -1. * crossthi

    crossthj = np.cross(rcj, -1 * fT)
    FFj = -1. * crossthj

    grain_list[i].an_a += FFi
    grain_list[j].an_a += FFj
