import network_ensemble as ne
import numpy as np
import newton_gmres as ng

def run_single_projective_step(network_ensemble, wait_interval, onmanifold_interval, collection_interval, projection_interval):
    network_ensemble.run(wait_interval)
    ncollect_intervals = onmanifold_interval/collection_interval
    coarse_evo = []
    coarse_evo.append(network_ensemble.restrict())
    for i in range(ncollect_intervals):
        network_ensemble.run(collection_interval)
        coarse_evo.append(network_ensemble.restrict())
    project(coarse_evo, projection_interval)

def project(coarse_evo, projection_interval):
    
