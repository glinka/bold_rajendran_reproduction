import network_model as nm
import numpy as np
import matplotlib.pyplot as plt

class Network_Ensemble:
    'ensemble of networks'

    def __init__(self, nnetworks, n, p, r):
        self.ensemble = []
        self.nnetworks = nnetworks
        self.n = n
        self.p = p
        self.r = r
        for i in range(nnetworks):
            self.ensemble.append(nm.Network_Model(n, p, r))
        
    def init_ensemble_havelhakimi(self, degs_seq):
        for network in self.ensemble:
            network.init_graph_havelhakimi(np.copy(degs_seq))

    def init_ensemble(self):
        for network in self.ensemble:
            network.init_er_graph()

    def run(self, nsteps):
        for network in self.ensemble:
            network.run(nsteps)
        
    def get_avg_deg_pdf(self, comm):
        rank = comm.Get_rank()
        nprocs = comm.Get_size()

        ensemble_degs = np.zeros((self.nnetworks, self.n))
        for i in range(self.nnetworks):
            ensemble_degs[i, :] = np.copy(self.ensemble[i].degs)

        all_ensembles = comm.gather(ensemble_degs, root=0)
        if rank == 0:
            ensemble_degs = np.concatenate(all_ensembles)
            ensemble_degs = np.sort(ensemble_degs.flatten())
            unique_degs = np.unique(ensemble_degs)
            nunique_degs = unique_degs.shape[0]
            deg_counts = np.zeros(nunique_degs)
            current_deg_index = 0
            for i in range(ensemble_degs.shape[0]):
                if ensemble_degs[i] == unique_degs[current_deg_index]:
                    deg_counts[current_deg_index] = deg_counts[current_deg_index] + 1
                else:
                    current_deg_index = current_deg_index + 1
                    deg_counts[current_deg_index] = deg_counts[current_deg_index] + 1

            # # pad with zeros for interpolation
            # deg_pdf = np.zeros(nunique_degs+1)
            # deg_pdf[0] = 0
            deg_pdf = deg_counts/float(np.sum(deg_counts))
            return unique_degs, deg_pdf
        else:
            return True, True


    def get_avg_deg_cdf(self, nbins=100):
        ensemble_degs = np.zeros((self.nnetworks, self.n))
        for i in range(self.nnetworks):
            ensemble_degs[i, :] = np.copy(self.ensemble[i].degs)
        ensemble_degs = np.sort(ensemble_degs.flatten())
        unique_degs = np.unique(ensemble_degs)
        nunique_degs = unique_degs.shape[0]
        deg_counts = np.zeros(nunique_degs)
        current_deg_index = 0
        for i in range(ensemble_degs.shape[0]):
            if ensemble_degs[i] == unique_degs[current_deg_index]:
                deg_counts[current_deg_index] = deg_counts[current_deg_index] + 1
            else:
                current_deg_index = current_deg_index + 1
                deg_counts[current_deg_index] = deg_counts[current_deg_index] + 1

        # pad with zeros for interpolation
        deg_pdf = np.zeros(nunique_degs+1)
        deg_pdf[0] = 0
        deg_pdf[1:] = deg_counts/float(np.sum(deg_counts))
        temp = np.copy(unique_degs)
        unique_degs = np.zeros(nunique_degs+1)
        unique_degs[1:] = temp
        unique_degs[0] = unique_degs[1] - 1


        # # TESTING
        # tdeg_counts = np.zeros(nunique_degs)
        # for i in range(nunique_degs):
        #     tdeg_counts[i] = sum(ensemble_degs == unique_degs[i+1])
        # print sum(tdeg_counts - deg_counts)
        # # END TESTING

        percentiles = np.arange(0,nbins+1)/float(nbins)
        deg_cdf = np.zeros(nbins)
        current_percentile = 0
        j = 0
        FLOAT_ERROR = 1e-8
        for i in range(nbins):
            while percentiles[i] > current_percentile + FLOAT_ERROR:
                j = j + 1
                current_percentile = current_percentile + deg_pdf[j]
            # linear interpolation time
            deg_cdf[i] = self.interpolate_deg(percentiles[i], unique_degs[j-1], np.sum(deg_pdf[:j-1]), unique_degs[j], np.sum(deg_pdf[:j]))
        return [deg_cdf, percentiles]


    def get_percentile_degs(self, comm, nbins=100):

        rank = comm.Get_rank()
        nprocs = comm.Get_size()

        ensemble_degs = np.zeros((self.nnetworks, self.n))
        for i in range(self.nnetworks):
            ensemble_degs[i, :] = np.copy(self.ensemble[i].degs)

        all_ensembles = comm.gather(ensemble_degs, root=0)
        if rank == 0:
            ensemble_degs = np.concatenate(all_ensembles)
            ensemble_degs = np.sort(ensemble_degs.flatten())
            unique_degs = np.unique(ensemble_degs)
            nunique_degs = unique_degs.shape[0]
            deg_counts = np.zeros(nunique_degs)
            current_deg_index = 0
            for i in range(ensemble_degs.shape[0]):
                if ensemble_degs[i] == unique_degs[current_deg_index]:
                    deg_counts[current_deg_index] = deg_counts[current_deg_index] + 1
                else:
                    current_deg_index = current_deg_index + 1
                    deg_counts[current_deg_index] = deg_counts[current_deg_index] + 1

            deg_pdf = deg_counts/float(np.sum(deg_counts))

            percentile_degs = np.empty(nbins+1)
            percentile_degs[0] = unique_degs[0] - 1

            current_percentile = 0.;
            index = 0;

            FLOAT_ERROR = 1e-8
            for i, percentile in enumerate(1./nbins*np.arange(1, nbins+1)):
                while current_percentile < percentile - FLOAT_ERROR:
                    current_percentile += deg_pdf[index]
                    index += 1
                percentile_degs[i+1] = unique_degs[index-1]
        else:
            percentile_degs = None

        percentile_degs = comm.bcast(percentile_degs, root=0)
        return percentile_degs


    def interpolate_deg(self, cdf_current, deg1, cdf1, deg2, cdf2):
        m = (cdf2 - cdf1)/float(deg2- deg1)
        return (cdf_current - (cdf2 - m*deg2))/m
        
    def project(self, deg_cdf_markers, times, proj_interval, comm):
        if comm.Get_rank() == 0:
            ncdfs = deg_cdf_markers.shape[0]
            nbins = deg_cdf_markers.shape[1]
            A = np.ones((ncdfs, 2))
            new_deg_cdf_markers = np.empty(nbins)
            for i in range(nbins):
                A[:,0] = times
                linear_coeffs = np.linalg.lstsq(A, deg_cdf_markers[:,i])[0]
                new_deg_cdf_markers[i] = linear_coeffs[0]*(times[-1] + proj_interval) + linear_coeffs[1]
                # # visual testing
                # fig = plt.figure()
                # ax = fig.add_subplot(111)
                # ax.scatter(times, deg_cdf_markers[:,i], lw=0, c='b')
                # ax.scatter(times[-1] + proj_interval, new_deg_cdf_markers[i], c='r')
                # ax.plot((times[0], times[-1] + proj_interval), linear_coeffs[0]*(np.array((times[0], times[-1] + proj_interval))) + linear_coeffs[1], c='g')
                # ax.set_xlim((times[0]-1, times[-1]+proj_interval+1))
                # plt.show()


            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # really shouldn't need to sort, nor does it really make sense to
            # but it makes things run
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            new_deg_cdf_markers = np.sort(new_deg_cdf_markers)

            # # visual testing
            # fig = plt.figure()
            # ax = fig.add_subplot(111)
            # for i in range(ncdfs):
            #     ax.plot(np.arange(nbins), deg_cdf_markers[i], color='b')
            # ax.plot(np.arange(nbins), new_deg_cdf_markers, color='r')
            # plt.show()
            

            test_network = nm.Network_Model(self.n, self.p, self.r)
            bin_width = 1./(nbins - 1)
            new_deg_percentiles = np.sort(np.random.uniform(size=self.n))
            new_degs = None
            graphical = False
            while not graphical:
                new_deg_percentiles = np.sort(np.random.uniform(size=self.n))
                new_degs = np.rint(new_deg_cdf_markers[np.rint(new_deg_percentiles/bin_width).astype(int)]).astype(int)
                graphical = test_network.init_graph_havelhakimi(new_degs[::-1])
        else:
            new_degs = None

        new_degs = comm.bcast(new_degs, root=0)
        return new_degs
