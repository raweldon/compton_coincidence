#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
import os
import glob

def single():
    ''' plot coincidence spectra
    '''
    cwd = os.getcwd()
    directory = cwd + '/plots/' 
    filenames = ('coinc_data_10ns.npy', 'coinc_data_20ns.npy')#, 'coinc_data_50ns.npy', 'coinc_data_100ns.npy') 

    for i, filename in enumerate(filenames):
        print filename

        ql_vs_tof = np.load(directory + filename)
        print ql_vs_tof
        ql = ql_vs_tof[0]
        print max(ql)
        tof = ql_vs_tof[1]

        '''Plots'''      
        plt.figure()
        plt.plot(tof, ql, 'o', alpha=0.5)
        plt.xlabel(r'$\Delta$t (ns)',fontsize=18)
        plt.ylabel('Light Output (MeVee)',fontsize=18)
        plt.title(os.path.basename(filename),fontsize=20)
        #plt.title(filename,fontsize=20)
        #plt.ylim(-0.2,6.0)
        #plt.savefig(save_dir + '/b_det'+str(det)+'.png', dpi=400)  

        plt.figure()
        plt.hist(ql, bins=50)
        plt.title(os.path.basename(filename),fontsize=20)
        
    plt.show()

def with_time():
    ''' plot edge spectra
    '''
    cwd = os.getcwd()
    directory = cwd + '/plots/' 
    filename = 'cs_hists_coinc_t2.npy'

    ql_vs_tof = np.load(directory + filename)
    ql = ql_vs_tof[0]

    indices = np.linspace(0, len(ql), 24)
    plt.figure()
    color = cm.viridis(np.linspace(0, 1, len(indices)))
    for idx, val in enumerate( indices):
        val = int(val)
        if idx == 0:
            ql_short = ql[0:val]
        else:
            ql_short = ql[int(indices[idx-1]):val]
        plt.hist(ql_short, bins=1000, histtype='step', color=color[idx], density=True)
        plt.title(os.path.basename(filename),fontsize=20)
        plt.xlim(0, 20000)
    plt.show()

if __name__ == '__main__':
    #single()
    with_time()
