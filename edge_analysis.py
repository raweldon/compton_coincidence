#!/usr/bin/env python
''' updated version of simultaneous_fit.py 
    fits uncalibrated ADC spectra to simulation
'''
import numpy as np
import matplotlib.pyplot as plt
import sys
import copy
from scipy.interpolate import interp1d
import lmfit
import os.path
import time

# directories
numpy_dir = '/media/radians/proraid/2018.1.tunl/beam_expt/out_numpy/'
sim_dir = '/home/radians/raweldon/tunl.2018.1_analysis/plastic_analysis/final_analysis/lo_calibration/simulated_spectra/'
prefix = 'beam_expt_'
run_database = 'run_database.txt'
start_time = time.time()

    
def get_meas_data(det_no, f_in, save_dir):   
    if det_no == 0:
        data = np.load(save_dir + f_in, allow_pickle=True)[0]
    else:
        data = np.load(save_dir + f_in, allow_pickle=True)[1]

    max_data = 80000    
    bin_width = max_data/500.
    meas_hist, meas_bin_edges = np.histogram(data, bins=np.arange(0, max_data + bin_width, bin_width))
    meas_bin_centers = (meas_bin_edges[:-1] + meas_bin_edges[1:])/2 
 
    meas_data = np.array((meas_bin_centers, meas_hist))

    #plt.figure()
    #plt.plot(meas_bin_centers, meas_hist)
    #plt.show()
    return meas_data    

def fit_range(low_bound, high_bound, data_in, name):
    '''select data for fit'''
    print '\nbounds:',name,low_bound,high_bound
    fit_range = (low_bound, high_bound)
    low_index = np.searchsorted(data_in[0], fit_range[0])
    high_index = np.searchsorted(data_in[0], fit_range[1])
    data_out = np.array((np.zeros(high_index-low_index), np.zeros(high_index-low_index)))
    data_out[0], data_out[1] = data_in[0][low_index:high_index], data_in[1][low_index:high_index]
    return data_out

def gaussian(x, mu, sigma, A):
    return A/np.sqrt(2.*np.pi*sigma**2) * np.exp(-(x - mu)**2/(2.*sigma**2))

def gaussian_smear(x, y, alpha, beta, gamma, c1, c2):
    '''smears bin centers (x) and bin contents (y) by a gaussian'''
    y_update = np.zeros(len(y))
    for i in xrange(len(x)):
        x_val = x[i]
        smear_val = np.sqrt((alpha*x_val)**2. + beta**2.*x_val + gamma**2.) 
        gaus_vals = gaussian(x, x_val, smear_val, 1)
        y_update[i] = np.dot(y, gaus_vals)
    return y_update
    
def shift_and_scale((alpha, beta, gamma, shift, spread, c1, c2, y_scale), bin_data, lin_scaling):
    '''
    shifts and scales the guassian
    bin_data is an (x,y) numpy array of bin centers and values
    shift allows for a horizontal shift in the bin centers (e.g. bin centers (0.5, 1, 1.5) --> (1, 1.5, 2) )
    spread allows for horizontal scaling in the bin centers (e.g. bin centers (0.5, 1, 1.5) --> (1, 2, 3) )
    '''
    x, y = bin_data
    y = gaussian_smear(x, y, alpha, beta, gamma, c1, c2)
    if lin_scaling:
        y *= c1
    else:
        y *= c1*x**c2
    x = copy.deepcopy(x)
    x = x*spread + np.ones(len(x))*shift
    return x, y

def lin_scaling(shift, spread, bin_data):
    ''' scales uncalibrated measured data '''
    x, y = bin_data
    x = copy.deepcopy(x)
    x = (x - np.ones(len(x))*shift)/spread
    return x, y
    
def calc_chisq(sim_x, sim_y, data_x, data_y):
    # bin centers in simulation are different from measurement -> use linear interpolation to compute SSE
    interp = interp1d(sim_x, sim_y, bounds_error=False, fill_value=0)
    sim_vals = interp(data_x)
    res = data_y - sim_vals
    chisq = [x**2/sim_vals[i]**2 for i, x in enumerate(res)]  
    return np.sqrt(chisq)

def minimize(fit_params, *args):
    '''
    Minimization based on length of input arguements 
    sim_data is an (x, y) of the simulated bin centers and bin contents
    meas_data is an (x, y) of the measured bin centers and contents
    '''
    pars = fit_params.valuesdict()
    alpha_1 = pars['alpha_1']
    beta_1 = pars['beta_1']
    gamma_1 = pars['gamma_1']
    c1_1 = pars['c1_1']
    c2_1 = pars['c2_1']
    shift = pars['shift']
    spread = pars['spread']

    y_scale_1 = pars['y_scale_1']
    sim_data, meas_data = args[0]
    sim_x_1, sim_y_1 = shift_and_scale((alpha_1, beta_1, gamma_1, shift, spread, c1_1, c2_1, y_scale_1), sim_data, args[1])
    data_x_1, data_y_1 = meas_data
    chisq_1 = calc_chisq(sim_x_1, sim_y_1, data_x_1, data_y_1)     

    return chisq_1
  
def spectra_fit(fit_params, *args, **kwargs):
    print '\nperforming minimization'

    fit_kws={'nan_policy': 'omit'}
    sim_data_1, meas_data_full_1, meas_data_1 = args[0]
    lin_scaling = args[1]
    print '    single sprectrum fit'
    res = lmfit.minimize(minimize, fit_params, method='leastsq', args=((sim_data_1, meas_data_1), lin_scaling), **fit_kws)

    if kwargs['print_info']:    
        print '\n',res.message
        print lmfit.fit_report(res)  
 
    if kwargs['show_plots']:
        plot_fitted_spectra( res.params['shift'].value, res.params['spread'].value, lin_scaling, 
                                ( (res.params['alpha_1'].value,), (res.params['beta_1'].value,), (res.params['gamma_1'].value,),
                                (res.params['c1_1'].value,), (res.params['c2_1'].value,),
                                (res.params['y_scale_1'].value,), (sim_data_1,), (meas_data_full_1,), (meas_data_1,) ) )

    # get shift term
    shift_term = res.params['shift'].value
    spread_term = res.params['spread'].value
    
    return shift_term, spread_term

def plot_fitted_spectra(shift, spread, lin_scale, args):
    for index, (alpha, beta, gamma, c1, c2, y_scale, sim_data, meas_data_full, meas_data) in enumerate(zip(args[0], args[1], 
                                                                                                           args[2], args[3], 
                                                                                                           args[4], args[5],
                                                                                                           args[6], args[7], args[8])):
        # update sim data
        sim_new = shift_and_scale((alpha, beta, gamma, shift, spread, c1, c2, y_scale), sim_data, lin_scale)
       
        # plot measured and fitted simulated data
        plt.figure()
        plt.plot(meas_data_full[0], meas_data_full[1], linestyle='None', marker='x', markersize=5, alpha=0.8, label='measured')
        plt.plot(sim_new[0], sim_new[1], '--', label='fit')
        plt.xlabel('ADC units', fontsize=20)
        plt.ylabel('Counts', fontsize=20)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlim(9000, 14500)
        plt.ylim(0, 30000)
        plt.legend(fontsize=18) 
        plt.tight_layout()
         
        # update measured data and plot
        meas_scaled = lin_scaling(shift, spread, meas_data_full)
        sim_with_res = shift_and_scale((alpha, beta, gamma, 0, 1, c1, c2, y_scale), sim_data, lin_scale)

        plt.figure()            
        plt.plot(meas_scaled[0], meas_scaled[1], linestyle='None', marker='x', markersize=5, alpha=0.5, label='measured')
        plt.plot(sim_with_res[0], sim_with_res[1], '--', label='sim with res')
        plt.plot(sim_data[0], sim_data[1], alpha=0.3, label='sim data')      
 
        # plot edge locations
        if index == 0:
            a = 0.3
            y_min, y_max = plt.ylim()
            y_range = np.arange(0, 5000 + 100, 100.)
            plt.plot([0.478]*len(y_range), y_range, 'k--', alpha=a)
            plt.text(0.478, y_max - y_max/15, 'cs edge')
            #plt.plot([0.343]*len(y_range), y_range, 'k--', alpha=a)
            #plt.text(0.343, y_max, 'na edge')
            #plt.plot([1.075]*len(y_range), y_range, 'k--', alpha=a)
            #plt.text(1.075, y_max - y_max/10,'na edge')
            #plt.plot([0.975]*len(y_range), y_range, 'k--', alpha=a)
            #plt.text(0.975, y_max, 'co edge')
            #plt.plot([1.133]*len(y_range), y_range, 'k--', alpha=a)
            #plt.text(1.133, y_max - y_max/5, 'co edge')

        plt.ylim(1, 80000)
        plt.xlim(0, 0.9)
        plt.xlabel('ql (MeVee)')
        plt.ylabel('Counts')
        plt.yscale('log')
        plt.legend()

def cal_interp(shift_terms, spread_terms, start_file, stop_file, run_names, database_dir):
    ''' Interpolates between calibrations to get calibration for each measurement.
        Pass start_file (beginning calibration) and stop_file (end_calibration).
        returns numpy array with groomed numpy file name and calibration factor.
    '''
    print run_names
    run = get_run_names([run_names], database_dir, run_database) # all files to parse
    groomed_arrays = []
    print '\nFiles between start and stop:'
    for i,r in enumerate(run):
        print r
        cal_names = ('cs','na','co')
        if any(cal_name in r for cal_name in cal_names):
               continue
        groomed_arrays.extend(get_numpy_arr(database_dir, run_database, r, numpy_dir, prefix, True))

    num_files = float(len(groomed_arrays))
    print '\nshift_terms, spread_terms:', shift_terms, spread_terms

    if shift_terms[0] == shift_terms[1]: # accounts from end of 11MeV run with no end calibration
        #print [shift_terms[0]]*len(groomed_arrays)
        return [shift_terms[0]]*len(groomed_arrays), [spread_terms[0]]*len(groomed_arrays), groomed_arrays

    else:
        step_shift = (max(shift_terms) - min(shift_terms))/(num_files - 1) 
        interp_shift = np.arange(min(shift_terms), max(shift_terms) + step_shift/2, step_shift)

        # account for positive or negative change in shift parameter
        if shift_terms[0] > shift_terms[1]:
            interp_shift = interp_shift[::-1]

        step_spread = (max(spread_terms) - min(spread_terms))/(num_files - 1) 
        interp_spread = np.arange(min(spread_terms), max(spread_terms) + step_spread/2, step_spread)

        if spread_terms[0] > spread_terms[1]:
            interp_spread = interp_spread[::-1]

        print 'interp shift[0], shift[-1]:', interp_shift[0], interp_shift[-1] 
        print 'interp spread[0], spread[-1]:', interp_spread[0], interp_spread[-1]
        return interp_shift, interp_spread, groomed_arrays


def single(fin, det_no, spread, min_range, lin_scaling):
    '''
    Use for individual spectrum fits to get rough idea of initial values and limits for simultaneous fitting
    Note: shift term is only accurate for simulataneous fit, do not use here
    '''

    iso = 'cs'
    cwd = os.getcwd()

    # load sim data
    if det_no == 0:
        ranges = [min_range, 14500] # stilbene
        sim_data = np.load(cwd + '/plots/cs_spec.npy')
        print '\nloading stilbene stimulated spectrum'
    if det_no == 1:
        ranges = [7000, 22000]
        sim_data = np.load(cwd + '/plots/ej309.npy')
        print '\nloading ej309 simulated spectrum'
    sim_data = [sim_data[0][20:], sim_data[1][20:]]
    print 'sim data loaded'

    meas_data_full = get_meas_data(det_no, fin, cwd + '/plots/')
    print 'meas data loaded'
    meas_data = fit_range(ranges[0], ranges[1], meas_data_full, iso)  #0.35,0.7
    data = [sim_data, meas_data_full, meas_data]            

    # lmfit (curve fit wrapper, default: leastsq Levenberg-Marquardt)
    # Only fit beta (Kornilov does this)
    fit_params = lmfit.Parameters()
    fit_params.add('alpha_1', value=0.0, min=0., max=20., vary=False)
    fit_params.add('beta_1', value=0.04, min=0.0035, max=0.1, vary=True)
    fit_params.add('gamma_1', value=0.0, min=0., max=20, vary=False)
    fit_params.add('shift', value=-1000, vary=True)
    fit_params.add('spread', value=spread, min=0, max=50000, vary=True) # 4mev 180000  
    if lin_scaling:
        fit_params.add('c1_1', value=0.01, vary=True)
        fit_params.add('c2_1', value=0., vary=False)
    else:     
        fit_params.add('c1_1', value=0.01, vary=True)
        fit_params.add('c2_1', value=-1.0, vary=True)
    fit_params.add('y_scale_1', value=0.01)
        
    e0, c = spectra_fit(fit_params, data, lin_scaling, print_info=True, show_plots=True)
    return e0, c


if __name__ == '__main__':
    ''' 
    '''
    fin = 'cs_hists_coinc_t2_2.npy'
    det_no = 0  # 0 stilbene, 1 ej309
    spread = 28000 # initial guess for spread
    min_range = 9000
    e0, c = single(fin, det_no, spread, min_range, lin_scaling=True)

    # print 477 keV estimate
    val = c*(0.477 + e0/c)
    print '\n477 keV ADC value =', round(val, 3)
    
    print ("--- %s seconds ---" % (time.time() - start_time))
    plt.show()
