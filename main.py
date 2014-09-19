import network_ensemble as ne
import numpy as np
import newton_gmres as ng
import sys
import time

def completion_bar(current, total, elapsed_time=None):
    perc = int((100.0*current)/total)
    bar = '\r['
    for i in range(perc):
        bar = bar + '|'
    for i in range(100-perc):
        bar = bar + ' '
    bar = bar + '] '
    if elapsed_time is not None:
        percf = (100.0*current)/total
        remaining_seconds = int(elapsed_time/(percf/100.0)) - int(elapsed_time)
        if remaining_seconds > 120:
            remaining_mins = int(100.0*remaining_seconds/60.0)/100.0
            bar = bar + str(remaining_mins) + ' m remaining'
        else:
            bar = bar + str(remaining_seconds) + 's remaining'
    print bar,
    sys.stdout.flush()

def main(nnetworks=20, nsteps=20000, n=100, p=0.5, r=0.9, project=True, proj_interval=1000, nintervals=50):

    ensemble = ne.Network_Ensemble(nnetworks, n, p, r)
    ensemble.init_ensemble()
    
    nbins = 50
    folder = None
    if project:
        folder = './py_data_project/'
    else:
        folder = './py_data_noproject/'
    files = {}
    files['deg_cdf_vals'] = open(folder + 'deg_probs.csv', 'w')
    files['degs'] = open(folder + 'degs.csv', 'w')
    files['times'] = open(folder + 'times.csv', 'w')
    for key in files.keys():
        files[key].write('filekey=' + key + '\n')

    coarse_vars = np.empty((nmicro_intervals, nbins))
    start = time.clock()

    if project:
        nmicrosteps = 1000
        micro_interval = 200
        nmicro_intervals = int(nmicrosteps/micro_interval)
        nintervals = int(nsteps/(nmicrosteps + proj_interval))
        interval = nmicrosteps + proj_interval
        times = np.empty(nmicro_intervals)
        print '-----------------------------------------'
        print '************** cpi enabled **************'
        print 'ensemble size:', nnetworks
        print 'n:', n
        print 'p:', p
        print 'r:', r
        print 'bin count:', nbins
        print 'total steps:', nsteps
        print 'micro interval:', micro_interval
        print 'projection interval:', proj_interval
        print 'total saves:', nintervals*nmicro_intervals
        print 'saving into: ', folder
        print '-----------------------------------------'
        for i in range(nintervals):
            for j in range(nmicro_intervals):
                ensemble.run(micro_interval)
                avg_deg_cdf = ensemble.get_avg_deg_cdf(nbins=nbins)
                coarse_vars[j,:] = avg_deg_cdf[0]
                times[j] = i*interval + (j+1)*micro_interval

                # save things
                avg_deg_cdf[0].shape = (1, avg_deg_cdf[0].shape[0])
                avg_deg_cdf[1].shape = (1, avg_deg_cdf[1].shape[0])
                np.savetxt(files['degs'], avg_deg_cdf[0], delimiter=',')
                np.savetxt(files['deg_cdf_vals'], avg_deg_cdf[1], delimiter=',')
                files['times'].write(str(interval*i + micro_interval*(j+1)) + '\n')

            new_deg_dist = ensemble.project(np.array(coarse_vars[-3:,:]), times[-3:], proj_interval)
            ensemble.init_ensemble_havelhakimi(new_deg_dist)
            completion_bar(i+1, nintervals, time.clock() - start)

    else:
        # save 100 avg deg cdfs
        interval = int(nsteps/nintervals)
        print '------------------------------------------'
        print '************** cpi disabled **************'
        print 'ensemble size:', nnetworks
        print 'n:', n
        print 'p:', p
        print 'r:', r
        print 'bin count:', nbins
        print 'total steps:', nsteps
        print 'save interval:', interval
        print 'total saves:', nintervals
        print 'saving into: ', folder
        print '------------------------------------------'
        for i in range(nintervals):
            ensemble.run(interval)
            avg_deg_cdf = ensemble.get_avg_deg_cdf(nbins=nbins)

            # save things
            avg_deg_cdf[0].shape = (1, avg_deg_cdf[0].shape[0])
            avg_deg_cdf[1].shape = (1, avg_deg_cdf[1].shape[0])
            np.savetxt(files['degs'], avg_deg_cdf[0], delimiter=',')
            np.savetxt(files['deg_cdf_vals'], avg_deg_cdf[1], delimiter=',')
            files['times'].write(str(interval*(i+1)) + '\n')
            completion_bar(i+1, nintervals, time.clock() - start)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', nargs=2, type=int, default=100)
    parser.add_argument('--nnetworks', type=int, default=20)
    parser.add_argument('--p', type=float, default=0.5)
    parser.add_argument('--r', type=float, default=0.9)
    parser.add_argument('--nintervals', type=int, default=50)
    parser.add_argument('--projectstep', type=int, default=1000)
    parser.add_argument('--nsteps', type=int, default=20000)
    parser.add_argument('--noproject', action='store_true', default=False)
    # not useful, but comforting
    parser.add_argument('--project', action='store_true', default=True)
    args = parser.parse_args()
    main(args.nnetworks, args.nsteps, args.n, args.p, args.r, project=not args.noproject, proj_interval=args.projectstep, nintervals=args.nintervals)
