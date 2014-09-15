import numpy as np

class Network_Model:
    'instance of the network model described in paper'

    def __init__(self, n=100, p=0.5, r=0.9):
        self.n = n
        self.p = p
        self.r = r
        self.A = np.zeros((self.n,self.n))
        self.degs = np.zeros(self.n)

    def init_er_graph(self):
        v = 1
        w = -1
        while v <= self.n:
            rv = np.random.uniform()
            w = int(w + 1 + np.floor(np.log(1 - rv)/np.log(1 - self.p)))
            while w >= v and v <= self.n:
                w = w - v
                v = v + 1
            if v < self.n:
                self.A[v][w] = 1
                self.A[w][v] = 1
                self.degs[v] = self.degs[v] + 1
                self.degs[w] = self.degs[w] + 1

        # BEGIN TESTING
        print 'initialization disparity:', np.sum(self.degs)/(self.n*(self.n-1)) - self.p
        # END TESTING
        
    def run(self, nsteps):
        v1_rns = (self.n*np.random.uniform(size=nsteps)).astype(int)
        v2_rns = (self.n*np.random.uniform(size=nsteps)).astype(int)
        edgeremoval_rns = np.random.uniform(size=nsteps)
        for i in range(nsteps):
            self.step(v1_rns[i], v2_rns[i], edgeremoval_rns[i])

    def step(self, v1, v2, edgeremoval_rn):
        if v1 != v2:
            if self.A[v1][v2] == 0:
                self.A[v1][v2] = 1
                self.A[v2][v1] = 1
                self.degs[v1] = self.degs[v1] + 1
                self.degs[v2] = self.degs[v2] + 1

            if edgeremoval_rn < self.r:
                removed_deg = int(np.sum(self.degs)*np.random.uniform())
                degcount = 0
                i = 0
                while degcount < removed_deg:
                    degcount = degcount + self.degs[i]
                    i = i + 1
                i = i - 1
                degcount = degcount - self.degs[i]
                j = 0
                while degcount < removed_deg:
                    degcount = degcount + self.A[i][j]
                    j = j + 1
                j = j - 1
                self.A[i][j] = 0
                self.A[j][i] = 0
                self.degs[i] = self.degs[i] - 1
                self.degs[j] = self.degs[j] - 1
