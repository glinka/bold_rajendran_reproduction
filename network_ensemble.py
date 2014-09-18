import network_model as nm
import numpy as np

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
        
    def get_avg_deg_cdf(self, nbins=10):
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

        percentiles = np.arange(1,nbins+1)/float(nbins)
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

    def interpolate_deg(self, cdf_current, deg1, cdf1, deg2, cdf2):
        m = (cdf2 - cdf1)/float(deg2- deg1)
        return (cdf_current - (cdf2 - m*deg2))/m
        
    def project(self, deg_cdf_markers, times, proj_interval):
        import matplotlib.pyplot as plt
        ncdfs = deg_cdf_markers.shape[0]
        nbins = deg_cdf_markers.shape[1]
        A = np.ones((ncdfs, 2))
        new_deg_cdf_markers = np.empty(nbins)
        for i in range(nbins):
            A[:,0] = times
            linear_coeffs = np.linalg.lstsq(A, deg_cdf_markers[:,i])[0]
            new_deg_cdf_markers[i] = linear_coeffs[0]*(times[-1] + proj_interval) + linear_coeffs[1]
            # TESTING
            # fig = plt.figure()
            # ax = fig.add_subplot(111)
            # ax.scatter(times, deg_cdf_markers[:,i], lw=0, c='b')
            # ax.scatter(times[-1] + proj_interval, new_deg_cdf_markers[i], c='r')
            # ax.plot((times[0], times[-1] + proj_interval), linear_coeffs[0]*(np.array((times[0], times[-1] + proj_interval))) + linear_coeffs[1], c='g')
            # ax.set_xlim((times[0]-1, times[-1]+proj_interval+1))
            # plt.show()
            # TESTING
        # pad front with start of cdf
        temp = np.copy(new_deg_cdf_markers)
        new_deg_cdf_markers = np.empty(nbins + 1)
        new_deg_cdf_markers[1:] = temp
        # kind of poor way of dealing with 
        new_deg_cdf_markers[0] = np.floor(temp[1])

        # create new, graphical degree sequence
        graphical = False
        new_degs = np.empty(self.n)
        test_network = nm.Network_Model(self.n, self.p, self.r)
        attempt_count = 0
        cdf_markers = np.arange(1,nbins+1)/float(nbins)
        test = None
        while not graphical:
            probs = np.random.uniform(size=self.n)
            for i in range(self.n):
                j = 0
                while probs[i] > cdf_markers[j]:
                    j += 1
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                # could also linearly interpolate between markers
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                new_degs[i] = int(new_deg_cdf_markers[j-1] + 0.5)
            if np.sum(new_degs) % 2 != 0:
                new_degs[-1] += 1
            # new_degs_test = test_network.init_graph_havelhakimi_test(np.sort(new_degs)[::-1].astype(int))
            # graphical = new_degs_test[0]
            graphical = test_network.init_graph_havelhakimi_test(np.sort(new_degs)[::-1].astype(int))
            # print attempt_count # np.sort(new_degs)[::-1].astype(int)
            attempt_count += 1
        print 'attempt count:', attempt_count, 'with a total of', np.sum(new_degs)/2, 'new edges'
        return new_degs
        # return new_degs_test[1]
