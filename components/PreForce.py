import numpy as np
import math
def preforce_grains(i, dt, grain_list):
    # global grain_list

    # v(t+dt/2)
    grain_list[i].v = 0.5 * grain_list[i].a * dt + grain_list[i].v
    grain_list[i].an_v = 0.5 * grain_list[i].an_a * dt + grain_list[i].an_v

    #x(t + dt)
    grain_list[i].pos = grain_list[i].pos + grain_list[i].v * dt
    grain_list[i].phi = grain_list[i].phi + grain_list[i].an_v * dt

    #reset accelerations before the calculation of forces
    grain_list[i].F = np.array([0., 0., 0.])
    grain_list[i].M = np.array([0., 0., 0.])
