import network_ensemble as ne
import numpy as np
import newton_gmres as ng

def main(nnetworks=20, nsteps=20000, n=100, p=0.5, r=0.9):

    ensemble = ne.Network_Ensemble(nnetworks, n, p, r)
    ensemble.init_ensemble()
    
    # save 100 avg deg dists
    nintervals = 50
    interval = int(nsteps/nintervals)
    deg_probs_out = open('./py_data/deg_probs.csv', 'w')
    degs_out = open('./py_data/degs.csv', 'w')
    times_out = open('./py_data/times.csv', 'w')
    for i in range(nintervals):
        print i/float(nintervals)
        ensemble.run(interval)
        avg_deg_dist = ensemble.get_avg_deg_dist(bins=10)

        avg_deg_dist[0].shape = (1, avg_deg_dist[0].shape[0])
        avg_deg_dist[1].shape = (1, avg_deg_dist[1].shape[0])
        np.savetxt(deg_probs_out, avg_deg_dist[0], delimiter=',')
        # output centers of bins, not edges
        bin_centers = (avg_deg_dist[1][0,:-1] + avg_deg_dist[1][0,1:])/2.0
        bin_centers.shape = (1, bin_centers.shape[0])
        np.savetxt(degs_out, bin_centers, delimiter=',')
        times_out.write(str(interval*(i+1)) + '\n')

    # necessary?
    deg_probs_out.close()
    degs_out.close()
    times_out.close()

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
