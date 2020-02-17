from coinc_analysis import fit_coinc_data
from edge_analysis import single
import os
import numpy as np
import time
import matplotlib.pyplot as plt

def main():

    # compton coinc
    bins = 300
    min_ql = 3000
    cwd = os.getcwd()
    directory = cwd + '/plots/' 
    filename = 'coinc_t2_4_10ns.npy'
    coinc_data = np.load(directory + filename)
    mu, sigma, counts = fit_coinc_data(coinc_data, min_ql, bins, gauss=True, show_plot=True)

    # compton edge fit
    fin = 'cs_hists_coinc_t2_3.npy'
    det_no = 0  # 0 stilbene, 1 ej309
    spread = 28000 # initial guess for spread
    min_range = 8000
    e0_p, c_p = single(fin, det_no, spread, min_range, lin_scaling=False)
    e0_l, c_l = single(fin, det_no, spread, min_range, lin_scaling=True)

    return mu, sigma, counts, e0_p, c_p, e0_l, c_l

if __name__ == '__main__':
    mu, sigma, counts, e0_p, c_p, e0_l, c_l = main()
    print '\n------------------------------------------'
    print '               Full results '
    print '------------------------------------------'
    print '\nCompton coincidence results:'
    print '{:^8s} {:>8s} {:>8s} {:^16s}'.format('mu', 'sigma', 'uncert', 'rel_uncert')
    print '{:^8.2f} {:>8.2f} {:>8.2f} {:>8.2f}%'.format(mu, sigma, sigma/np.sqrt(counts), sigma/mu/np.sqrt(counts)*100)

    print '\nCompton edge results:'
    val_l = c_l*(0.477 + e0_l/c_l)
    val_p = c_p*(0.477 + e0_p/c_p)
    print '{:^10s} {:^15s} {:^14s} {:^12s}'.format('linear', 'lin-coinc diff', 'powerlaw', 'power-coinc diff')
    print '{:^8.2f} {:^15.2f}% {:^14.2f} {:^12.2f}%'.format(val_l, abs(val_l - mu)/mu*100, val_p, abs(val_p - mu)/mu*100)
    plt.show()

