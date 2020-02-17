import numpy as np
import os
import matplotlib.pyplot as plt
import polimi_parse

cwd = os.getcwd()
sim_dir = cwd + '/plots/'
input_file = 'cs_spec_polimi.log'

collision = polimi_parse.fromfile(sim_dir + input_file)
histories = np.unique(collision.history)
energy_deposition=[]

for history in histories:
    #print history
    unique_hist_array = collision[np.where(collision.history == history)]
    delta_e = np.sum(unique_hist_array.deltaEnergy)
    energy_deposition.append(delta_e)
energy_deposition = np.asarray(energy_deposition)

max_data = 1.5
bin_width = max_data/500. # changing this doesn't seem to make a difference (1/31/19)
#make bin centers in the same x location for any spectrum
sim_hist, sim_bin_edges = np.histogram(energy_deposition,
                                    bins=np.arange(0, max_data + bin_width, bin_width))
sim_bin_centers = (sim_bin_edges[:-1] + sim_bin_edges[1:])/2
sim_data = np.array((sim_bin_centers, sim_hist))

plt.figure()
plt.plot(sim_bin_centers, sim_hist)
plt.show()

np.save(sim_dir + 'cs_spec', sim_data) # save to npy for quick loading
print 'spectra saved to ' + sim_dir + 'ej309'
