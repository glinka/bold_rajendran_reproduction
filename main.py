import network_ensemble as ne
import numpy as np
import newton_gmres as ng
import sys
import time
from mpi4py import MPI
import matplotlib.pyplot as plt

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

def main(nnetworks=100, nsteps=20000, n=100, p=0.5, r=0.9, project=True, proj_interval=1000):

    # init mpi
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nprocs = comm.Get_size()
    nnetworks_per_proc = nnetworks/nprocs

    ensemble = ne.Network_Ensemble(nnetworks_per_proc, n, p, r)
    ensemble.init_ensemble()

    nbins = 100
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

    start = time.clock()

    if project:
        micro_interval = 100
        nmicro_intervals = 10
        nmicrosteps = micro_interval*nmicro_intervals
        interval = nmicrosteps + proj_interval
        nintervals = 2
        coarse_vars = np.empty((nmicro_intervals, nbins+1))
        times = np.empty(nmicro_intervals)
        if rank == 0:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            print '-----------------------------------------'
            print '************** cpi enabled **************'
            print 'ensemble size:', nnetworks_per_proc*nprocs
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
                percentile_degs = ensemble.get_percentile_degs(comm)

                avg_deg_pdf, degs = ensemble.get_avg_deg_pdf(comm)
                if rank == 0:
                    ax.plot(avg_deg_pdf, degs, color='b')

                coarse_vars[j] = percentile_degs
                times[j] = i*interval + (j+1)*micro_interval

            new_deg_dist = ensemble.project(np.array(coarse_vars[-5:]), times[-5:], proj_interval, comm)
            ensemble.init_ensemble_havelhakimi(new_deg_dist)
            if rank == 0:
                completion_bar(i+1, nintervals, time.clock() - start)

        if rank == 0:
            plt.show()

    else:


        interval = 1000
        nintervals = nsteps/interval

        avg_deg_pdf, degs = ensemble.get_avg_deg_pdf(comm)

        if rank == 0:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            # capture initial conditions
            ax.plot(avg_deg_pdf, degs, color='b')

            print '------------------------------------------'
            print '************** cpi disabled **************'
            print 'ensemble size:', nnetworks_per_proc*nprocs
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
            # avg_deg_cdf = ensemble.get_avg_deg_cdf(nbins=nbins)

            # # save things
            # avg_deg_cdf[0].shape = (1, avg_deg_cdf[0].shape[0])
            # avg_deg_cdf[1].shape = (1, avg_deg_cdf[1].shape[0])
            # np.savetxt(files['degs'], avg_deg_cdf[0], delimiter=',')
            # np.savetxt(files['deg_cdf_vals'], avg_deg_cdf[1], delimiter=',')
            # files['times'].write(str(interval*(i+1)) + '\n')
            # completion_bar(i+1, nintervals, time.clock() - start)

            avg_deg_pdf, degs = ensemble.get_avg_deg_pdf(comm)

            if rank == 0:
                ax.plot(avg_deg_pdf, degs, color='b')

            # # save things
            # avg_deg_pdf[0].shape = (1, avg_deg_cdf[0].shape[0])
            # avg_deg_cdf[1].shape = (1, avg_deg_cdf[1].shape[0])
            # np.savetxt(files['degs'], avg_deg_cdf[0], delimiter=',')
            # np.savetxt(files['deg_cdf_vals'], avg_deg_cdf[1], delimiter=',')
            # files['times'].write(str(interval*(i+1)) + '\n')
            # completion_bar(i+1, nintervals, time.clock() - start)

        if rank == 0:
            plt.show()

if __name__=='__main__':
    # import argparse
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--n', nargs=2, type=int, default=100)
    # parser.add_argument('--nnetworks', type=int, default=20)
    # parser.add_argument('--p', type=float, default=0.5)
    # parser.add_argument('--r', type=float, default=0.9)
    # parser.add_argument('--nintervals', type=int, default=50)
    # parser.add_argument('--projectstep', type=int, default=1000)
    # parser.add_argument('--nsteps', type=int, default=20000)
    # parser.add_argument('--noproject', action='store_true', default=False)
    # # not useful, but comforting
    # parser.add_argument('--project', action='store_true', default=True)
    # args = parser.parse_args()
    # main(args.nnetworks, args.nsteps, args.n, args.p, args.r, project=not args.noproject, proj_interval=args.projectstep, nintervals=args.nintervals)
    main()
