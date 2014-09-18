import argparse
import numpy as np
import matplotlib.pyplot as plt

import util_fns as uf

def deg_cdf_evo(data, ax=None, cmap='jet'):
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    plot = False
    if ax is None:
        plot = True
        fig = plt.figure()
        ax = fig.add_subplot(111)
    times = data['times']
    ntimes = times.shape[0]
    colornorm = colors.Normalize(vmin=times[0], vmax=times[-1])
    colormap = cm.ScalarMappable(norm=colornorm, cmap=cmap)
    for i in range(ntimes):
        ax.plot(data['degs'][i,:], data['deg_cdf_vals'][i,:], c=colormap.to_rgba(1.0*times[i]))
    if plot:
        plt.show()

def compare_deg_cdf(project_data, noproject_data):
    # own method in util_fns faster? need to get matching index set of both arrays
    matching_project_indices = np.where(np.in1d(project_data['times'], noproject_data['times'], assume_unique=True))
    matching_noproject_indices = np.where(np.in1d(noproject_data['times'], project_data['times'], assume_unique=True))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for key in project_data.keys():
        project_data[key] = project_data[key][matching_project_indices]
    for key in noproject_data.keys():
        noproject_data[key] = noproject_data[key][matching_noproject_indices]
    deg_cdf_evo(project_data, ax=ax, cmap='Blues')
    deg_cdf_evo(noproject_data, ax=ax, cmap='Reds')
    plt.show()
    
                


def main(args):
    import os
    noproject_data = {}
    project_data = {}
    for filename in args.input_files:
        params = uf.get_header_data(open(os.path.realpath(filename), 'r').readline())
        filekey = params['filekey']
        if 'noproject' in filename:
            noproject_data[filekey], garbage = uf.get_data(filename, header_rows=1)
        elif '_project' in filename:
            project_data[filekey], garbage = uf.get_data(filename, header_rows=1)
        else:
            print 'input filename not recognized:', filename
    if args.compare_deg_cdf:
        compare_deg_cdf(project_data, noproject_data)
    if args.deg_cdf:
        if project_data:
            deg_cdf_evo(project_data)
        if noproject_data:
            deg_cdf_evo(noproject_data)
            
            

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input_files', nargs='+')
    parser.add_argument('--compare-deg-cdf', action='store_true', default=False)
    parser.add_argument('--deg-cdf', action='store_true', default=False)
    args = parser.parse_args()
    main(args)


