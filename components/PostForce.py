import numpy as np
import math

def postforce_grains(gravity,i, dt, grain_list,endForce):
    # global endForce

    #no gravity in this test case

    #compute acceleration a=f/m
    grain_list[i].a = grain_list[i].F/grain_list[i].mass()
    grain_list[i].an_a = grain_list[i].M/ grain_list[i].inertia()

    #velocity v(t + dt)
    grain_list[i].v = 0.5 * grain_list[i].a * dt + grain_list[i].v
    grain_list[i].an_v = 0.5 * grain_list[i].an_a * dt + grain_list[i].an_v
