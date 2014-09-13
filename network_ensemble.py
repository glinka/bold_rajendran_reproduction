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
        
    def init_ensemble(self):
        for network in self.ensemble:
            network.init_er_graph()

    def run(self, nsteps):
        for network in self.ensemble:
            network.run(nsteps)
        
    def get_avg_deg_dist(self):
        ensemble_degs = np.zeros((self.nnetworks, self.n))
        for i in range(self.nnetworks):
            ensemble_degs[i, :] = np.copy(self.ensemble[i].degs)
        maxdeg = np.amax(ensemble_degs)
        mindeg = np.amin(ensemble_degs)
        print maxdeg, mindeg
        return np.histogram(ensemble_degs, bins=maxdeg-mindeg+1, range=(mindeg, maxdeg), density=True)[0]
        
        
