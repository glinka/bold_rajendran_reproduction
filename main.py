import network_ensemble as ne
import numpy as np

def main(nnetworks=20, nsteps=30000, n=100, p=0.5, r=0.9):

    ensemble = ne.Network_Ensemble(nnetworks, n, p, r)
    ensemble.init_ensemble()
    
    # save 100 avg deg dists
    nintervals = 100
    interval = int(nsteps/nintervals)
    f = open('./py_data/deg_dists.csv', 'w')
    for i in range(nintervals):
        print i/float(nintervals)
        ensemble.run(interval)
        avg_deg_dist = ensemble.get_avg_deg_dist()
        avg_deg_dist.shape = (1, avg_deg_dist.shape[0])
        np.savetxt(f, avg_deg_dist, delimiter=',')

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', nargs=2, type=int, default=100)
    parser.add_argument('--nnetworks', type=int, default=20)
    parser.add_argument('--p', type=float, default=0.5)
    parser.add_argument('--r', type=float, default=0.9)
    parser.add_argument('--nsteps', type=int, default=30000)
    args = parser.parse_args()
    print args
    main(args.nnetworks, args.nsteps, args.n, args.p, args.r)
