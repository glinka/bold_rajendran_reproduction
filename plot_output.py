import argparse
import numpy as np
import matplotlib.pyplot as plt

def deg_dist_evo(times, degs, deg_probs):
    import matplotlib.colorbar as cb
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ntimes = times.shape[0]
    print degs.shape
    print deg_probs.shape
    for i in range(ntimes):
        ax.plot(degs[i,:], deg_probs[i,:], c=(float(i)/ntimes, 0.5, 0.5))
    plt.show()

def main():
    deg_dist_evo(np.genfromtxt('./py_data/times.csv', delimiter=','), np.genfromtxt('./py_data/degs.csv', delimiter=','), np.genfromtxt('./py_data/deg_probs.csv', delimiter=','))

if __name__=='__main__':
    main()
