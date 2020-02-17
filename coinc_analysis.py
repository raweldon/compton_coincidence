import numpy as np 
import os
import lmfit
import matplotlib.pyplot as plt
from scipy.special import erfc,erf

def gaussian(x, a, mu, sig):
    return a*np.exp( -(x - mu)**2.0/(2.0*sig**2.0) )

def brubaker(x, mu, sig, m, b):
    return ((m*x + b)/2. * (1. - erf((x - mu)/(sig*np.sqrt(2.)))) -
            m*sig/np.sqrt(2.*np.pi)*np.exp(-(x - mu)**2./(2.*sig**2.)))

def gaus_exp_convo_line(x, a, mu, sig, tau, m, b):
    erfc_arg = (sig/tau - (x-mu)/sig)/np.sqrt(2.)
    func = a*sig/tau * np.sqrt(np.pi/2.) * np.exp(-(0.5*(sig/tau)**2 - (x-mu)/tau)) * (1-erf(erfc_arg))
    return func 

def gauss_line(x, a, mu, sig, m, b):
    ''' gaussian plus brubakers convolution'''
    gaus = a*np.exp( -(x - mu)**2.0/(2.0*sig**2.0) )
    #line = m*x + b
    error = ((m*x + b)/2. * (1. - erf((x - mu)/(sig*np.sqrt(2.)))) -
            m*sig/np.sqrt(2.*np.pi)*np.exp(-(x - mu)**2./(2.*sig**2.)))
    return gaus + error

def fit_coinc_data(coinc_data, min_ql, bins, gauss, show_plot):
    ql, tof = coinc_data
    ql = ql[np.where(ql > min_ql)]
    ql_hist, bins = np.histogram(ql, bins=bins, range=(ql.min(), 19000))
    bin_centers = (bins[:-1] + bins[1:])/2

    if gauss:
        # fit histogram 
        gmodel = lmfit.Model(gauss_line)
        params = gmodel.make_params(a=10, mu=13500, sig=400, m=1, b=0)
        params['mu'].max = 13700
        params['mu'].min = 13200
        #params['b'].max = -1000
        #params['b'].vary = False
        #params['m'].max = 0
        params['sig'].max = 800
        params['sig'].min = 0
        res = gmodel.fit(ql_hist, params=params, x=bin_centers, nan_policy='omit')#, method='nelder')
        print '\n', lmfit.fit_report(res)

        # get number of counts under gaussian (estimate of uncert on mean)
        mu = res.params['mu'].value
        sigma = res.params['sig'].value
        counts = len(ql[np.where((ql < mu + 3*sigma ) & (ql > mu - 3*sigma))])

        if show_plot:
            plt.figure(100)
            #plt.plot(bin_centers, ql_hist, linestyle='none', marker='x', alpha=0.8)
            plt.plot(bin_centers, ql_hist, alpha=0.9, label='measured')
            plt.plot(bin_centers, gauss_line(bin_centers, res.params['a'], res.params['mu'], res.params['sig'], res.params['m'], res.params['b']), 
                     '--', linewidth=2.5, zorder=100, label='fit')
            #plt.plot([mu - 3*sigma]*100, np.linspace(0, 60, 100), '--k')
            #plt.plot([mu + 3*sigma]*100, np.linspace(0, 60, 100), '--k')
            plt.plot(bin_centers, brubaker(bin_centers, res.params['mu'], res.params['sig'], res.params['m'], res.params['b']), 
                     '-.', linewidth=2.5, zorder=100, label='background fit')
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            plt.ylabel('Counts', fontsize=20)
            plt.xlabel('ADC units', fontsize=20)
            plt.legend(fontsize=18)
            plt.tight_layout()

        return res.params['mu'].value, res.params['sig'].value, counts

    else:
        # fit histogram 
        gmodel = lmfit.Model(gaus_exp_convo_line)
        params = gmodel.make_params(a=50, mu=13500, sig=400, tau=1, m=0, b=0)
        params['mu'].max = 13700
        params['mu'].min = 13200
        #params['b'].min = 1
        #params['m'].max = 0
        params['sig'].max = 800
        params['sig'].min = 0
        res = gmodel.fit(ql_hist, params=params, x=bin_centers, nan_policy='omit')#, method='nelder')
        print '\n', lmfit.fit_report(res)

        # get number of counts under gaussian (estimate of uncert on mean)
        mu = res.params['mu'].value
        sigma = res.params['sig'].value
        counts = len(ql[np.where((ql < mu + 3*sigma ) & (ql > mu - 3*sigma))])

        if show_plot:
            plt.figure(100)
            plt.plot(bin_centers, ql_hist)
            plt.plot(bin_centers, gauss_line(bin_centers, res.params['a'], res.params['mu'], res.params['sig'], res.params['m'], res.params['b']))
            plt.plot([mu - 3*sigma]*100, np.linspace(0, 20, 100), '--r')
            plt.plot([mu + 3*sigma]*100, np.linspace(0, 20, 100), '--r')
            
        
        return res.params['mu'].value + res.params['tau'].value, np.sqrt(res.params['sig'].value**2 + res.params['tau']**2), counts


if __name__ == '__main__':
    ''' set gauss=True for gauss fit, flase for gaussian exponential convolution
    '''

    cwd = os.getcwd()
    directory = cwd + '/plots/' 
    filename = 'coinc_t2_4_10ns.npy'

    bins = 300
    ql_min = 4000
    coinc_data = np.load(directory + filename)
    mu, sigma, counts = fit_coinc_data(coinc_data, ql_min, bins, gauss=True, show_plot=True)
    print '\n mu =', round(mu, 3), '\n sigma =', round(sigma, 3), '\n counts =', counts
    plt.show()