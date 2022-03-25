import numpy as np
import math
import os
from components.ParticlesClass import Particles
from components.BondClass import Bond
from components.IdentifyBonds import identify_cohesive_bonds
from components.BondKernel import bonded_grains_interaction
from components.DryKernel import dryParticle
from components.PreForce import preforce_grains
from components.PostForce import postforce_grains
from components.VTK import write_vtk_file
from components.Delta import delta_calculation
from components.Plot import plotDelta



bond_dic={}
grain_list = []
E = 2e11
ForceGlobal = 100000
en = 0.2
ks_to_kn = 1.
nus_to_nun = 0.2
mu = 0.3
gravity = np.array([0.,0.,0.])
R = 0.2
rhoS = 2600.
dt = 1e-5
endForce = np.array([0.,0.,0.])
timesteps = 20000
vtkInterval=500


def particlesGenerator():
    positions = []
    positions = []
    positions.append(np.array([0., 0., 0.]))
    for i in range(1, 11):
        positions.append(np.array([2 * R * i - 0 * i, 0., 0.]))

    for i in range(len(positions)):
        prt = Particles(int(i), R, positions[i], rhoS)
        grain_list.append(prt)

def terminalInfo(tstep,vtkInterval ):
    if (tstep % vtkInterval == 0):
        print("timestep: " + str(tstep))
        print("endForce: " + str(endForce))

def preForceGrains():
    for index in range(len(grain_list)):
        preforce_grains(index, dt, grain_list)

def particlesInteraction():
    for index1 in range(len(grain_list)):
        for index2 in range(index1 + 1, len(grain_list)):
            if ((grain_list[index1].id, grain_list[index2].id) in bond_dic.keys()):
                bonded_grains_interaction(index1, index2, dt, grain_list, bond_dic)
def postForceGrains():
    for index in range(len(grain_list)):
        postforce_grains(gravity, index, dt, grain_list, endForce)

def fixParticle(i=0):
    grain_list[i].v = np.array([0., 0., 0.])
    grain_list[i].an_v = np.array([0., 0., 0.])
    grain_list[i].pos = np.array([0., 0., 0.])
    grain_list[i].phi = np.array([0., 0., 0.])

def vtkWriter(tstep):
    if (tstep % vtkInterval == 0):
        write_vtk_file(grain_list, tstep, r"VTKs\particles",100)
def forceIncrement(tstep, fIncrement, i=10):
    global endForce
    if (tstep <= timesteps // 2):
        endForce = endForce + fIncrement
    grain_list[i].F = grain_list[i].F +  endForce


def showForce(i, step, interval):
    if(step%interval == 0):
        print(f"particle {i} force is: {grain_list[i].F}")

def plotData(tstep):

    if (tstep % 4000 == 0):
        [delta_ys, delta_bs, L0] = delta_calculation(endForce, grain_list, bond_dic)
        plotDelta(delta_ys, delta_bs, L0, tstep, grain_list, bond_dic)

def main():
    global endForce
    global ForceGlobal
    Force = ForceGlobal

    particlesGenerator()
    identify_cohesive_bonds(grain_list, E, R, rhoS, en,bond_dic)
    test = bond_dic
    fIncrement = [0.,Force/(timesteps/2.),0.]

    for tstep in range(timesteps):

        # terminalInfo(tstep, vtkInterval)
        showForce(10, tstep, 100)
        preForceGrains()
        fixParticle(0)
        forceIncrement(tstep, fIncrement, 10)
        particlesInteraction()
        postForceGrains()
        vtkWriter(tstep)
        if(tstep%1000 == 0):
         delta_y, delta_b,L0 = delta_calculation(endForce, grain_list, bond_dic)
         plotDelta(delta_y, L0, tstep,grain_list)



if __name__ == "__main__":
    # execute only if run as a script

    main()
    print("finished!")
