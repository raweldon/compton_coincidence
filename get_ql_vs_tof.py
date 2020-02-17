#!/usr/bin/env python

'''
    Updated version from daq2 for compton coincidence analysis
    stil - det0
    ej309 - det 1
'''

import numpy as np
import time
import sys
import directories
from get_numpy_arrays import get_run_names, get_numpy_arr
import os

start_time = time.time()

def get_coinc_counts(names, database_dir, run_database, numpy_dir, window):

    runs = get_run_names(names, database_dir, run_database)
    for index,r in enumerate(runs):
        print '\n', r

        #cal_files = ['cs', 'na', 'co']
        #if any(cal in r for cal in cal_files):
        #    print 'Skipping calibration file:', r
        #    continue

        # parse data for ql and tof
        print numpy_dir + r
        data_files = get_numpy_arr(database_dir, run_database, r, numpy_dir, 'comp_coinc_', groomed=True)
        
        ql_dets, del_t_dets, qs_dets, backing_idx_dets = [], [], [], []
        for data_index, datum in enumerate(data_files):
            print data_index, datum
            f_in = np.load(numpy_dir + datum)       
            data = f_in['data']             
            scatter_index = np.where(data['det_no']==0)[0]
            
            det_index = np.where(data['det_no']==1)
            det_data = data[det_index[0]]
            backing = det_index[0]
                
            # search for scatter to backing hits with delta t < 100ns
            j=0
            i=0
            delta_t, ql, qs, backing_idx = [], [], [], []
            while j<(len(scatter_index)-1) and i<(len(backing)-1):
                if (scatter_index[j]<backing[i]):
                    coincidence = data['trig_offset'][backing[i]] + data['trig_time'][backing[i]] - \
                                  data['trig_offset'][scatter_index[j]] - data['trig_time'][scatter_index[j]] 
                    if coincidence < window:
                        delta_t.append(coincidence) 
                        ql.append(data['ql'][scatter_index[j]])
                    j+=1
                else:
                    i+=1
            
            ql_dets.extend(ql)
            del_t_dets.extend(delta_t)
     
            np.save(save_dir + 'coinc_data_' + str(window) + 'ns', (ql_dets, del_t_dets))

if __name__ == '__main__':
    ''' 
    '''
    window = 10 # ns, coincidence window
    numpy_dir, root_dir, save_dir, database_dir, run_database = directories.compton_dirs()
    names = ['coinc_t1']

    print 'saving data to', save_dir
    get_coinc_counts(names, database_dir, run_database, numpy_dir, window)
        
print("--- %s seconds ---" % (time.time() - start_time))


