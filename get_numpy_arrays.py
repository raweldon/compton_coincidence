import os
import sys

''' Contents:
    get_run_names - get file names for run_database.txt
    get_start_stop - get start and stop timestamps from run_database.txt
    get_numpy_arr - get numpy arrays by retrieving timestamps 
    get_root_files - get root files by retrieving timestamps
'''

def get_run_names(names,database_dir,run_database):
    ''' Get run name from run_database.txt
        Pass an array of file names to look for
    '''
    run=[]
    with open(database_dir + run_database, 'r') as infile:
        for line in infile:
            if 'comment' in line:
                break
            if any(name in line for name in names):
                run.append(line[1:-2])
        if not run:
            print 'Did not find files in run_database.txt'
            sys.exit()
    return run 

def get_start_stop(database_dir,filename,run):
    #get start and stop time stamps
    f = open(database_dir+filename)
    flag = True
    while flag == True:
        line = f.readline()
        tmp = line.strip()
        if not tmp or tmp[0] != '@':
            continue
        tmp = tmp[1:-1]
        if tmp == run:
            while tmp[0] != '}':
                line = f.readline()
                tmp = line.strip()
                tmp1 = line[0:line.find('=')]
                if tmp1 == 'start':
                    start = line[line.find('=')+1:-1]
                    #print 'start = ', start
                elif tmp1 == 'stop':
                    stop = line[line.find('=')+1:-1]
                    #print 'stop = ', stop
                else:
                    flag = False
    return start,stop

def get_numpy_arr(database_dir,filename,run,numpy_dir,prefix,groomed):
    '''
    get_numpy_arr:
    returns a list of strings of numpy array files for parsing.  Database_dir is a string 
    with the file dir, filename is run_database.txt (or name of file containing 
    the run info), run is the name of the files to be parsed from run_database.txt,
    numpy_dir is the location of the numpy arrays.
    '''
    start, stop = get_start_stop(database_dir,filename,run)
           
    #get files for parsing
    numpy_arrays = os.listdir(numpy_dir)        
    numpy_arrays = sorted(numpy_arrays, key=lambda item: (int(item.partition('_')[0])
                          if item[0].isdigit() else float('inf'), item))
    flag = True
    i = 0
    names=[]
    while flag == True:
        if groomed == True:
            # check if raw_cleaned files
            if 'raw.npz' in numpy_arrays[0].split('_'):
                if numpy_arrays[i] == prefix+start+'_groomed_raw.npz':
                    while numpy_arrays[i] != prefix+stop+'_groomed_raw.npz':
                        if 'no_beam' in numpy_arrays[i]:
                            i = i + 1
                            continue
                        names.append(numpy_arrays[i])
                        i = i + 1
                    names.append(prefix+stop+'_groomed_raw.npz')
                    flag = False
                i = i + 1
            else:    
                if numpy_arrays[i] == prefix+start+'_groomed.npz':
                    while numpy_arrays[i] != prefix+stop+'_groomed.npz':
                        if 'no_beam' in numpy_arrays[i]:
                            i = i + 2
                            continue
                        names.append(numpy_arrays[i])
                        i = i + 2
                    names.append(prefix+stop+'_groomed.npz')
                    flag = False
                i = i + 1
        else:
            if numpy_arrays[i] == prefix+start+'.npz':
                while numpy_arrays[i] != prefix+stop+'.npz':
                    names.append(numpy_arrays[i])
                    i = i + 2
                names.append(prefix+stop+'.npz')
                flag = False
            i = i + 1

    return names

def get_numpy_arr_tof(database_dir,filename,run,numpy_dir,prefix):
    ''' same as get_numpy_arr but for 0 deg tof analysis
    '''
    start, stop = get_start_stop(database_dir,filename,run)
           
    #get files for parsing
    numpy_arrays = os.listdir(numpy_dir)        
    numpy_arrays = sorted(numpy_arrays, key=lambda item: (int(item.partition('_')[0])
                          if item[0].isdigit() else float('inf'), item))
    flag = True
    i = 0
    data=[]
    while flag == True:
        if numpy_arrays[i] == prefix+start+'_0.npz':
            while numpy_arrays[i] != prefix+stop+'_9.npz': # when tof files split into 10 separate files
                data.append(numpy_arrays[i])
                i = i + 1
            data.append(prefix+stop+'_9.npz')
            flag = False
        i = i + 1
    return data



def get_root_files(database_dir,filename,run,root_dir,prefix):
    '''
    same as get_numpy_arr but for root tree files
    '''
    start, stop = get_start_stop(database_dir,filename,run)            
    #get files for parsing
    tree = os.listdir(root_dir)        
    tree = sorted(tree, key=lambda item: (int(item.partition(' ')[0])
                          if item[0].isdigit() else float('inf'), item))
    flag = True
    i = 0
    data=[]
    while flag == True:
        if tree[i] == prefix+start+'.root':
            while tree[i] != prefix+stop+'.root':
                data.append(tree[i])
                i = i + 1
            data.append(prefix+stop+'.root')
            flag = False
        i = i + 1

    return data



