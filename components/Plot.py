import matplotlib.pyplot as plt
import numpy as np
import math
def plotDelta(delta_ys, L0, step,grain_list):
    # global grain_list
    # global bond_dic
    xs = [x.pos[0] for x in grain_list]
    ys = [y.pos[1] for y in grain_list]
    ys = np.array(ys)
    xs = [x/L0 for x in xs]
    fig = plt.figure(figsize=(5, 5))
    plt.xlabel('x/L0')
    plt.ylabel('displacement (mm)')
    plt.title("for timestep: " + str(step))
    plt.plot(xs, delta_ys*1000, 'o', label='delta_ys')
    plt.plot(xs, ys*1000, 'o',label ='simulation')
    plt.legend(loc='best')
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    plt.savefig("./Plots/"+str(step)+'.png')
    plt.close()
