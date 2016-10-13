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

        # # BEGIN TESTING
        # print 'initialization disparity:', np.sum(self.degs)/(self.n*(self.n-1)) - self.p
        # # END TESTING
        
    def init_graph_havelhakimi(self, deg_sequence):
        # assumes reverse-sorted deg sequence as input, DOES NOT SORT WITHIN FUNCTION
        if np.sum(deg_sequence) % 2 == 1:
            deg_sequence[0] += 1
        self.A = np.zeros((self.n,self.n))
        self.degs = np.zeros(self.n)
        for i in range(self.n):
            current_deg = deg_sequence[i]
            j = i+1
            while current_deg > 0 and j < self.n:
                if deg_sequence[j] > 0:
                    self.A[i][j] = 1
                    self.A[j][i] = 1
                    self.degs[i] = self.degs[i] + 1
                    self.degs[j] = self.degs[j] + 1
                    deg_sequence[j] -= 1
                    current_deg -= 1
                j += 1
            deg_sequence[i] = 0
        if sum(deg_sequence == 0) == self.n:
            return True
        else:
            return False
            
    @staticmethod
    def init_graph_havelhakimi_test(deg_sequence):
        # assumes reverse-sorted deg sequence as input, DOES NOT SORT WITHIN FUNCTION
        # adjusted_deg_sequence = np.copy(deg_sequence)
        n = deg_sequence.shape[0]
        THRESHOLD = 0
        ncut_degs = 0
        for i in range(n):
            current_deg = deg_sequence[i]
            # while current_deg + i > n - 1:
            if current_deg + i > n - 1:
                return False
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                # if the difference is less than the threshold, decrease
                # degree to make graphical
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                # this needs fixing, as we were previously returning some sort of strange
                # zeroed deg_sequence which somehow still gave somewhat reasonable results
                # if we want to use this approach, need to make copy of sequence and adjust
                # appropriately, then return the adjusted copy
                # for the time being, simply test graphical/not graphical
                # if current_deg + i - THRESHOLD > n - 1:
                #     # print i, current_deg + i - (n - 1)
                #     return [False]
                # else:
                #     ncut_degs += current_deg - (n-1-i)
                #     deg_sequence[i] = n-1-i
                #     adjusted_deg_sequence[i] = n-1-i
                #     deg_sequence[i:] = np.sort(deg_sequence[i:])
                #     adjusted_deg_sequence[i:] = np.sort(adjusted_deg_sequence[i:])
                #     current_deg = deg_sequence[i]
            for j in range(i + 1, current_deg + i + 1):
                deg_sequence[i] = deg_sequence[i] - 1
                deg_sequence[j] = deg_sequence[j] - 1
        if sum(deg_sequence == 0) == n:
            # print 'havelhakimi cut', ncut_degs, 'degs for a total of', np.sum(deg_sequence)/2, 'edges'
            return True #, deg_sequence]
        else:
            print 'havelhakimi deg sum discrepancy:', sum(deg_sequence == 0) - n
            return False
        

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
    
        
