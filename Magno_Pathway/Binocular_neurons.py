#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 14:52:49 2017

@author: marcorax93

This script it has been developed on spyder, it is devided in sometimes
indipendent subsections, marked with #%%, thus do not run the entire script as whole

Using spyder is advised, also remember that the plots make use of pyqtgraph, 
for its speed on intensive 3d plotting, for this reason remember to change
the backend from inline to automatic in (Tools/preferences/IPython console/Graphics)
 
"""

# General computation and plotting
from pyqtgraph.Qt import QtCore
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import numpy as np
import scipy.io as sio
import pickle
import time
from pathlib import Path
from brian2 import amp, ms, mV, pA, mA, nS, nA, pF, us, volt, second, prefs,\
    SpikeMonitor, StateMonitor, figure, plot, show, xlabel, ylabel,\
    seed, xlim, ylim, subplot, network_operation, TimedArray,\
    defaultclock, SpikeGeneratorGroup, asarray, pamp, set_device, device
# Teili lib
from teili.core.groups import Neurons, Connections
from teili import teiliNetwork
from teili.models.neuron_models import DPI 
from teili.models.synapse_models import DPISyn
# My libs
from Tools.Stimuli import rotating_bar
from Tools.Connectivity import stereo_centers, connect_spike_gen, gaussian_spike_gen


from CellsParameters.fast_inhibitory_neuron_param import fast_in_parameters
from CellsParameters.fast_inhibitory_synapse_param import fast_in_parameters_syn
from CellsParameters.slower_neuron import slower_parameters


# Settings for c++ code building, it helps to speed up the computation a lot in
# long simulations 
# Note that build on run is set to false, meaning that the network won't rebuild 
# when calling.run() (In this way is possible to change some parameters and run
# the network without losing time to rebuild)

prefs.codegen.target = "numpy"
set_device('cpp_standalone', directory='C++ Buids/Rotating_Artificial_Bar', build_on_run=False)
prefs.devices.cpp_standalone.extra_make_args_unix = ['-j 16'] #The total number of file being built in parallel at once (this value depend on the number of cores, schedulers and so on) tweak it to find the best value on your machine
prefs.devices.cpp_standalone.openmp_threads = 8 #The maximum number of threads the run will feed 

#The default clock of the simulation
defaultclock.dt = 1 * ms

#Important to find the Path to the Data folder
parent_folder=str(Path().resolve().parent)


#%% Net Parameters

# This values can be changed without recompiling the network (standalone params)
scal_weight_disparity_neu=300 # This is a scaling weight applied to all disparity sensitive neurons
scal_weight_inhibitory_neu=600 # This is a scaling weight applied to all inhibitory neurons
weight_inhibitory = -1320 # This is the weight for inhibitory synapses


input_size_x=40
#input_size_y=80# Recording
input_size_y=40 # Artificial rotating bar
num_orientations=8
num_input_neurons=input_size_x*input_size_y

# Connectivity params 
kernel_size=5 #linear dimension of the squared gaussian connectivity kernel
sigma=1. 
spacing=2 # the number of pixel separating the centers of each consecutive
          # kernel, the bigger the number the smaller, the resolution of the 
          # encoded disparity
          
depth_planes=5 # The number of depth planes encoded, the bigger, higher is 
               # the maximum disparity encoded by the network, but accordingly,
               # the net will have more neurons and synapses

neuron_centers, num_stereo_cells = stereo_centers(input_size_x, input_size_y, spacing,
                                          kernel_size, depth_planes)

#%% Load simple cells data from previous simulations
# In the case of the artificial stimuli I will load the same file and then 
# apply a disparity in the next section
## L
inputFile = parent_folder+'/Data/Magno_Population_Data/Rotating_Artificial_Bar.data'
#inputFile = parent_folder+'/Data/Magno_Population_Data/L_On_Moving_Bar.data'
fd = open(inputFile, 'rb')
pop_data_L = pickle.load(fd)
## R
inputFile = parent_folder+'/Data/Magno_Population_Data/Rotating_Artificial_Bar.data'
#inputFile = parent_folder+'/Data/Magno_Population_Data/R_On_Moving_Bar.data'
fd = open(inputFile, 'rb')
pop_data_R = pickle.load(fd)



#%% Artificial stimuli disparity
# In the case of an artificial stimuli, the offset or disparity is set here
# Useful to check network ability to encode disparity
# Be careful to check that the offset data will stay within network boundaries

xoff = 2
yoff = 3
for k in range(num_orientations):
    pop_data_R[str(k)] = pop_data_R[str(k)] + xoff + yoff*input_size_x


#%% Spike generators - Simple Neruon stimuli
# Set two spike generators that will play the data recorded during the esecution
# of Monocular_neurons.py
    
simple_neurons_L=[]
simple_neurons_R=[]
for k in range(num_orientations):
    simple_neurons_L.append(SpikeGeneratorGroup(num_input_neurons, pop_data_L[str(k)], pop_data_L[str(k)+'t']))
    simple_neurons_R.append(SpikeGeneratorGroup(num_input_neurons, pop_data_R[str(k)], pop_data_R[str(k)+'t']))
     

#%% Setting up Network
    
start_time = time.time()

start = pop_data_L['start']
end = pop_data_L['end']

Net = teiliNetwork()

# Neurons
disparity_neurons = []
inhibitory_neurons_L = []
inhibitory_neurons_R = []
# Below I also set up the connectivity kernel centers for each neuron, the 
# position is referred to the input neuron spatial position, so for a disparity
# neuron with (xl,yl) = (10,10) (xr,yr) = (15,10) will have two gaussian connectivity
# kernel centered in position (10,10) in the left simple neuron sub plane
# and in position (15,20) in the right simple neuron sub plane, and it will have 
# the highest response to a disparity of (15-10,10-10) = (5,0)
# Please note that each disparity neuron is selective to a given orientation
# since it is connected to subplanes with the same orientation

# Inhibitory neurons take the name from the excitatory sublayer to which are connected
# inhibitory_neurons_L are excited from left monocular neurons and inhibited 
# from right monocular neurons
# inhibitory_neurons_R are excited from right monocular neurons and inhibited 
# from left monocular neurons
for k in range(num_orientations):
    disparity_neurons.append(Neurons(N=num_stereo_cells, name="disparity_neurons_pop_"+str(k), equation_builder=DPI(num_inputs=4)))
#    disparity_neurons[k].set_params(slower_parameters)
    disparity_neurons[k].variables.add_array('xl', constant=True, scalar=False, size=num_stereo_cells)
    disparity_neurons[k].variables.add_array('yl', constant=True, scalar=False, size=num_stereo_cells)
    disparity_neurons[k].variables.add_array('xr', constant=True, scalar=False, size=num_stereo_cells)
    disparity_neurons[k].variables.add_array('yr', constant=True, scalar=False, size=num_stereo_cells)
    disparity_neurons[k].xl = neuron_centers[:,0]
    disparity_neurons[k].yl = neuron_centers[:,1]
    disparity_neurons[k].xr = neuron_centers[:,2]
    disparity_neurons[k].yr = neuron_centers[:,3]
    
    inhibitory_neurons_L.append(Neurons(N=num_stereo_cells, name="inhibitory_neurons_L_pop_"+str(k), equation_builder=DPI(num_inputs=4)))
#    inhibitory_neurons_L[k].set_params(fast_in_parameters)   
    inhibitory_neurons_L[k].variables.add_array('xl', constant=True, scalar=False, size=num_stereo_cells)
    inhibitory_neurons_L[k].variables.add_array('yl', constant=True, scalar=False, size=num_stereo_cells)
    inhibitory_neurons_L[k].variables.add_array('xr', constant=True, scalar=False, size=num_stereo_cells)
    inhibitory_neurons_L[k].variables.add_array('yr', constant=True, scalar=False, size=num_stereo_cells)
    inhibitory_neurons_L[k].xl = neuron_centers[:,0]
    inhibitory_neurons_L[k].yl = neuron_centers[:,1]
    inhibitory_neurons_L[k].xr = neuron_centers[:,2]
    inhibitory_neurons_L[k].yr = neuron_centers[:,3]    
    
    inhibitory_neurons_R.append(Neurons(N=num_stereo_cells, name="inhibitory_neurons_R_pop_"+str(k), equation_builder=DPI(num_inputs=4)))
#    inhibitory_neurons_R[k].set_params(fast_in_parameters) 
    inhibitory_neurons_R[k].variables.add_array('xl', constant=True, scalar=False, size=num_stereo_cells)
    inhibitory_neurons_R[k].variables.add_array('yl', constant=True, scalar=False, size=num_stereo_cells)
    inhibitory_neurons_R[k].variables.add_array('xr', constant=True, scalar=False, size=num_stereo_cells)
    inhibitory_neurons_R[k].variables.add_array('yr', constant=True, scalar=False, size=num_stereo_cells)
    inhibitory_neurons_R[k].xl = neuron_centers[:,0]
    inhibitory_neurons_R[k].yl = neuron_centers[:,1]
    inhibitory_neurons_R[k].xr = neuron_centers[:,2]
    inhibitory_neurons_R[k].yr = neuron_centers[:,3]


# Synapses
# Input sinapses connecting disparity and inhibitory neurons to the simple monocular
# neurons from the subplanes underneath (input to the network for brevity)
disparity_neu_input_synapses_L = []
disparity_neu_input_synapses_R = []
inibitory_neu_L_input_synapses_L = []
inibitory_neu_L_input_synapses_R = []
inibitory_neu_R_input_synapses_L = []
inibitory_neu_R_input_synapses_R = []

# Synapses connecting inhibitory neurons to dipsarity neurons
inibitory_neu_L_disparity_neu_syn = []
inibitory_neu_R_disparity_neu_syn = []


# At the moment the only way to load functions with addictional parameters
# it's through string sobstitution, I don't know if it's a tili or brian2 problem
synapses_params = {'input_size_x': input_size_x,
                   'input_size_y': input_size_y,
                   'kernel_size': kernel_size,
                   'sigma': sigma}


for k in range(num_orientations):
    
    disparity_neu_input_synapses_L.append(Connections(simple_neurons_L[k], disparity_neurons[k],
                     name="disparity_neu_input_synapses_L_pop_"+str(k), equation_builder=DPISyn()))
    disparity_neu_input_synapses_R.append(Connections(simple_neurons_R[k], disparity_neurons[k],
                     name="disparity_neu_input_synapses_R_pop_"+str(k), equation_builder=DPISyn()))
        
    inibitory_neu_L_input_synapses_L.append(Connections(simple_neurons_L[k], inhibitory_neurons_L[k],
                     name="inibitory_neu_L_input_synapses_L_pop_"+str(k), equation_builder=DPISyn()))
    inibitory_neu_L_input_synapses_R.append(Connections(simple_neurons_R[k], inhibitory_neurons_L[k],
                     name="inibitory_neu_L_input_synapses_R_pop_"+str(k), equation_builder=DPISyn()))
    inibitory_neu_R_input_synapses_L.append(Connections(simple_neurons_L[k], inhibitory_neurons_R[k],
                     name="inibitory_neu_R_input_synapses_L_pop_"+str(k), equation_builder=DPISyn()))
    inibitory_neu_R_input_synapses_R.append(Connections(simple_neurons_R[k], inhibitory_neurons_R[k],
                     name="inibitory_neu_R_input_synapses_R_pop_"+str(k), equation_builder=DPISyn()))

    inibitory_neu_L_disparity_neu_syn.append(Connections(inhibitory_neurons_L[k], disparity_neurons[k],
                     name="inibitory_neu_L_disparity_neu_syn_pop_"+str(k), equation_builder=DPISyn()))  
    inibitory_neu_R_disparity_neu_syn.append(Connections(inhibitory_neurons_R[k], disparity_neurons[k],
                     name="inibitory_neu_R_disparity_neu_syn_pop_"+str(k), equation_builder=DPISyn()))



    disparity_neu_input_synapses_L[k].connect('connect_spike_gen(i, int(xl_post), int(yl_post), {input_size_x}, {input_size_y}, {kernel_size})'.format(**synapses_params))
    disparity_neu_input_synapses_R[k].connect('connect_spike_gen(i, int(xr_post), int(yr_post), {input_size_x}, {input_size_y}, {kernel_size})'.format(**synapses_params))
    inibitory_neu_L_input_synapses_L[k].connect('connect_spike_gen(i, int(xl_post), int(yl_post), {input_size_x}, {input_size_y}, {kernel_size})'.format(**synapses_params))
    inibitory_neu_L_input_synapses_R[k].connect('connect_spike_gen(i, int(xr_post), int(yr_post), {input_size_x}, {input_size_y}, {kernel_size})'.format(**synapses_params))
    inibitory_neu_R_input_synapses_L[k].connect('connect_spike_gen(i, int(xl_post), int(yl_post), {input_size_x}, {input_size_y}, {kernel_size})'.format(**synapses_params))
    inibitory_neu_R_input_synapses_R[k].connect('connect_spike_gen(i, int(xr_post), int(yr_post), {input_size_x}, {input_size_y}, {kernel_size})'.format(**synapses_params))

    inibitory_neu_L_disparity_neu_syn[k].connect(condition='i==j')
    inibitory_neu_R_disparity_neu_syn[k].connect(condition='i==j')
#    inibitory_neu_L_disparity_neu_syn[k].set_params(fast_in_parameters_syn)
#    inibitory_neu_R_disparity_neu_syn[k].set_params(fast_in_parameters_syn)

    disparity_neu_input_synapses_L[k].namespace['gaussian_spike_gen'] = gaussian_spike_gen
    disparity_neu_input_synapses_R[k].namespace['gaussian_spike_gen'] = gaussian_spike_gen
    inibitory_neu_L_input_synapses_L[k].namespace['gaussian_spike_gen'] = gaussian_spike_gen
    inibitory_neu_L_input_synapses_R[k].namespace['gaussian_spike_gen'] = gaussian_spike_gen
    inibitory_neu_R_input_synapses_L[k].namespace['gaussian_spike_gen'] = gaussian_spike_gen
    inibitory_neu_R_input_synapses_R[k].namespace['gaussian_spike_gen'] = gaussian_spike_gen

        
    disparity_neu_input_synapses_L[k].add_state_variable('scal_weight_disparity_neu',shared = True, constant = True, changeInStandalone=True)
    disparity_neu_input_synapses_L[k].scal_weight_disparity_neu = scal_weight_disparity_neu
    disparity_neu_input_synapses_L[k].weight = 'scal_weight_disparity_neu*gaussian_spike_gen(i, int(xl_post), int(yl_post), {input_size_x}, {input_size_y}, {sigma})'.format(**synapses_params)

    disparity_neu_input_synapses_R[k].add_state_variable('scal_weight_disparity_neu',shared = True, constant = True, changeInStandalone=True)
    disparity_neu_input_synapses_R[k].scal_weight_disparity_neu = scal_weight_disparity_neu
    disparity_neu_input_synapses_R[k].weight = 'scal_weight_disparity_neu*gaussian_spike_gen(i, int(xr_post), int(yr_post), {input_size_x}, {input_size_y}, {sigma})'.format(**synapses_params)

    inibitory_neu_L_input_synapses_L[k].add_state_variable('scal_weight_inhibitory_neu',shared = True, constant = True, changeInStandalone=True)
    inibitory_neu_L_input_synapses_L[k].scal_weight_inhibitory_neu = scal_weight_inhibitory_neu
    inibitory_neu_L_input_synapses_L[k].weight = 'scal_weight_inhibitory_neu*gaussian_spike_gen(i, int(xl_post), int(yl_post), {input_size_x}, {input_size_y}, {sigma})'.format(**synapses_params)

    inibitory_neu_L_input_synapses_R[k].add_state_variable('scal_weight_inhibitory_neu',shared = True, constant = True, changeInStandalone=True)
    inibitory_neu_L_input_synapses_R[k].scal_weight_inhibitory_neu = -scal_weight_inhibitory_neu # This connectivity kernel is inhibitory
    inibitory_neu_L_input_synapses_R[k].weight = 'scal_weight_inhibitory_neu*gaussian_spike_gen(i, int(xr_post), int(yr_post), {input_size_x}, {input_size_y}, {sigma})'.format(**synapses_params)    

    inibitory_neu_R_input_synapses_L[k].add_state_variable('scal_weight_inhibitory_neu',shared = True, constant = True, changeInStandalone=True)
    inibitory_neu_R_input_synapses_L[k].scal_weight_inhibitory_neu = -scal_weight_inhibitory_neu # This connectivity kernel is inhibitory
    inibitory_neu_R_input_synapses_L[k].weight = 'scal_weight_inhibitory_neu*gaussian_spike_gen(i, int(xl_post), int(yl_post), {input_size_x}, {input_size_y}, {sigma})'.format(**synapses_params)

    inibitory_neu_R_input_synapses_R[k].add_state_variable('scal_weight_inhibitory_neu',shared = True, constant = True, changeInStandalone=True)
    inibitory_neu_R_input_synapses_R[k].scal_weight_inhibitory_neu = scal_weight_inhibitory_neu 
    inibitory_neu_R_input_synapses_R[k].weight = 'scal_weight_inhibitory_neu*gaussian_spike_gen(i, int(xr_post), int(yr_post), {input_size_x}, {input_size_y}, {sigma})'.format(**synapses_params) 
        
    inibitory_neu_L_disparity_neu_syn[k].add_state_variable('weight_inhibitory',shared = True, constant = True, changeInStandalone=True)    
    inibitory_neu_L_disparity_neu_syn[k].weight_inhibitory = weight_inhibitory
    inibitory_neu_L_disparity_neu_syn[k].weight = "weight_inhibitory"
    
    inibitory_neu_R_disparity_neu_syn[k].add_state_variable('weight_inhibitory',shared = True, constant = True, changeInStandalone=True)    
    inibitory_neu_R_disparity_neu_syn[k].weight_inhibitory = weight_inhibitory
    inibitory_neu_R_disparity_neu_syn[k].weight = "weight_inhibitory"

#Monitors

input_spikes_L = []
input_spikes_R = []
stereo_spikes = []

# Spiking activity
for k in range(num_orientations):       
    input_spikes_L.append(SpikeMonitor(simple_neurons_L[k], name='spikemon_simple_neurons_L_pop_'+str(k)))
    input_spikes_R.append(SpikeMonitor(simple_neurons_R[k], name='spikemon_simple_neurons_R_pop_'+str(k)))
    stereo_spikes.append(SpikeMonitor(disparity_neurons[k],  name='spikemon_disparity_neurons_pop_'+str(k)))



# This monitor records states of an entire "Corolla"
# This sub population is composed by all the neurons sensitive to different
# disparities but with the same left connectivity kernel.

# "Position of the Corolla"
x0 = 20
y0 = 20

left_centers_bool_array = (neuron_centers[:,0]==x0)*(neuron_centers[:,1]==y0)

ind = np.where(left_centers_bool_array==True)[0]

disparity_neurons_corolla = []
# Also ihibitory neurons are stereo, therefore they can be represented 
# with Corollas
inhibitory_neurons_L_corolla = []
inhibitory_neurons_R_corolla = []

for k in range(num_orientations):       
    disparity_neurons_corolla.append(StateMonitor(disparity_neurons[k], variables=["Imem"], record=ind, name='statemon_disparity_neurons_corolla_pop_'+str(k)))
    inhibitory_neurons_L_corolla.append(StateMonitor(inhibitory_neurons_L[k], variables=["Imem"], record=ind, name='statemon_inhibitory_neurons_L_corolla_pop_'+str(k)))
    inhibitory_neurons_R_corolla.append(StateMonitor(inhibitory_neurons_R[k], variables=["Imem"], record=ind, name='statemon_inhibitory_neurons_R_corolla_pop_'+str(k)))
        
    
Net.add(disparity_neurons, inhibitory_neurons_L, inhibitory_neurons_R, simple_neurons_L,
        simple_neurons_R, disparity_neu_input_synapses_L, disparity_neu_input_synapses_R,
        inibitory_neu_L_input_synapses_L, inibitory_neu_L_input_synapses_R,
        inibitory_neu_R_input_synapses_L, inibitory_neu_R_input_synapses_R,
        inibitory_neu_L_disparity_neu_syn, inibitory_neu_R_disparity_neu_syn,
        stereo_spikes, input_spikes_L, input_spikes_R, disparity_neurons_corolla,
        inhibitory_neurons_L_corolla, inhibitory_neurons_R_corolla)

Net.build()

print("Network complete building time :" + str(time.time() - start_time))


#%% Run without compilation
# Setting standalone params, if you want to test different parameters,
# ovverride this ones and than after you find the right ones change them at the 
# beginning of the file 

for k in range(num_orientations):
    Net.standalone_params['disparity_neu_input_synapses_L_pop_'+str(k)+'_scal_weight_disparity_neu'] = scal_weight_disparity_neu 
    Net.standalone_params['disparity_neu_input_synapses_R_pop_'+str(k)+'_scal_weight_disparity_neu'] = scal_weight_disparity_neu 
    Net.standalone_params['inibitory_neu_L_input_synapses_L_pop_'+str(k)+'_scal_weight_inhibitory_neu'] = scal_weight_disparity_neu 
    Net.standalone_params['inibitory_neu_L_input_synapses_R_pop_'+str(k)+'_scal_weight_inhibitory_neu'] = -scal_weight_disparity_neu 
    Net.standalone_params['inibitory_neu_R_input_synapses_L_pop_'+str(k)+'_scal_weight_inhibitory_neu'] = -scal_weight_disparity_neu 
    Net.standalone_params['inibitory_neu_R_input_synapses_R_pop_'+str(k)+'_scal_weight_inhibitory_neu'] = scal_weight_disparity_neu 
    Net.standalone_params['inibitory_neu_L_disparity_neu_syn_pop_'+str(k)+'_scal_weight_inhibitory_neu'] = weight_inhibitory 
    Net.standalone_params['inibitory_neu_R_disparity_neu_syn_pop_'+str(k)+'_scal_weight_inhibitory_neu'] = weight_inhibitory

Net.run(end-start)


# Output Managing
# the results are packed here, and it will be used for plotting and saving
disparity_neurons_activity = {}
for k in range(num_orientations):
    disparity_neurons_activity[str(k)] = neuron_centers[stereo_spikes[k].i[:]]
    disparity_neurons_activity[str(k)+'t'] = stereo_spikes[k].t[:]
disparity_neurons_activity['start'] = start
disparity_neurons_activity['end'] = end


#%% Visualization Frametiser - Moving Bar
# Load the stimuli used to produce the Monocular_neuron.py activity,
# In this way you will visualize the original stimuli in the 3D visualizer
stimuli_res_x = 240
stimuli_res_y = 180

window_size_x = input_size_x
window_size_y = input_size_y

inputFile = parent_folder+'/Data/Magno_Population_Data/L_Moving_Bar_input.data'
fd = open(inputFile, 'rb')
Left_Retina = pickle.load(fd)

inputFile = parent_folder+'/Data/Magno_Population_Data/R_Moving_Bar_input.data'
fd = open(inputFile, 'rb')
Right_Retina = pickle.load(fd)

indexarray_input_L=Left_Retina[0][0].astype(int)
indexarray_input_R=Right_Retina[0][0].astype(int)

timearray_input_L=Left_Retina[0][1]*ms
timearray_input_R=Right_Retina[0][1]*ms

#%% Visualization Frametiser - Rotating Bar
# Load the stimuli used to produce the Monocular_neuron.py activity,
# In this way you will visualize the original stimuli in the 3D visualizer
stimuli_res_x = 240
stimuli_res_y = 180

window_size_x = input_size_x
window_size_y = input_size_y
# This settings need to be the same used to produce the stimuli for the monocular
# population
length = 15
center_offx = 0
center_offy = 0
angle_step = 40
ts_offset = 10 * defaultclock.dt

left_bar = rotating_bar(length, stimuli_res_x, stimuli_res_y, center_offx,
                                 center_offy, angle_step, ts_offset)[1]

indexarray_input_L = (left_bar[0] + left_bar[1]*stimuli_res_x).astype(int)
indexarray_input_R = (left_bar[0] + left_bar[1]*stimuli_res_x + xoff + yoff*stimuli_res_x).astype(int)
timearray_input_L = left_bar[2]*ts_offset
timearray_input_R = left_bar[2]*ts_offset

#%% Frametizer 
# In order to display an animated output of the network decoded disparity
# I need to produce frames of events.
# It also computes a vector displaying the direction of the stimuli decoded 
# in the window. NB it will take a while! Have a second coffe


timeWindow = 10 * defaultclock.dt # The width in time used to build a frame out of single events

counter = 0 * defaultclock.dt

#I want to display each population with a different color
Popcolors = [[1, 0, 0], [89/255, 198/255, 118/255], [0, 0, 1], [247/255, 0, 255/255],
          [0, 1, 0], [1, 128/255, 0], [120/255, 120/255, 120/255], [0, 171/255, 1]]

InposL = []
InposR = []
IncolorsL = []
IncolorsR = []
OverallDisparty = []
Disparities = np.zeros((len(neuron_centers),2))
Disparities[:,0] = neuron_centers[:,2]-neuron_centers[:,0]
Disparities[:,1] = neuron_centers[:,3]-neuron_centers[:,1]
# All the position for the Corollas
corollas_centers = np.unique(neuron_centers[:,0:2],axis=0)
SingleDisparities = []
while ((end-start)>counter):
    #Left 
    timetmp = (timearray_input_L>=counter)*(timearray_input_L<counter+timeWindow)
    z = np.zeros(len(timearray_input_L[timetmp==1]))
    tmpindex=np.unravel_index(indexarray_input_L[timetmp==1], (stimuli_res_y, stimuli_res_x))    
    ## the x coordinates are already right, the y coordinate have to be mirrored. The coordinate references are different from the frame to the plot 
    InposL.append(np.transpose(np.asarray([tmpindex[:][1], -tmpindex[:][0] + stimuli_res_y, z], dtype=np.float32 )))
    
    color = np.empty((len(z), 4), dtype=np.float32)
    color[:, 3] = 0.9
    color[:, 0] = np.clip(0, 0, 1)
    color[:, 1] = np.clip(1, 0, 1)
    color[:, 2] = np.clip(0, 0, 1)
    IncolorsL.append(color)
    #Right
    timetmp = (timearray_input_R>=counter)*(timearray_input_R<counter+timeWindow)
    z = np.zeros(len(timearray_input_R[timetmp==1]))
    tmpindex=np.unravel_index(indexarray_input_R[timetmp==1], (stimuli_res_y, stimuli_res_x))    
    ## the x coordinates are already right, the y coordinate have to be mirrored. The coordinate references are different from the frame to the plot 
    InposR.append(np.transpose(np.asarray([tmpindex[:][1], -tmpindex[:][0] + stimuli_res_y, z], dtype=np.float32 )))
    
    color = np.empty((len(z), 4), dtype=np.float32)
    color[:, 3] = 0.9
    color[:, 0] = np.clip(1, 0, 1)
    color[:, 1] = np.clip(0, 0, 1)
    color[:, 2] = np.clip(0, 0, 1)
    IncolorsR.append(color)
    tmpDisparity = np.asarray([0,0])
    tmpSingleDisparities = [] #all the disparities for all the cells in the populations in a single timewindow
    for k in range(num_orientations):
        timetmp = (stereo_spikes[k].t>counter)*(stereo_spikes[k].t<=counter+timeWindow)
        tmpindex=stereo_spikes[k].i[timetmp==1]
        tmpCenters = disparity_neurons_activity[str(k)][timetmp==1]
        tmpDisparities = np.zeros((len(corollas_centers),2))
        if len(tmpindex) != 0 :
            tmpDisparity = tmpDisparity + (sum(Disparities[tmpindex,:])/len(tmpindex))*np.asarray([1,-1])#The y coordinate must be reversed before plortting
            count=0
            for kk in corollas_centers:
                #Matricial Way to find all the kk in tmpCenters
                kkCentered_firing = (np.abs((tmpCenters[:,0:2]-kk))[:,0]+np.abs((tmpCenters[:,0:2]-kk))[:,1])==0
                if sum(kkCentered_firing) != 0:
                    tmpDisparities[count] = np.asarray([tmpCenters[kkCentered_firing][:,2]-tmpCenters[kkCentered_firing][:,0],
                                                tmpCenters[kkCentered_firing][:,3]-tmpCenters[kkCentered_firing][:,1]]).sum(axis=1)/(sum(kkCentered_firing)+len(tmpindex)/10)
                count = count + 1               #The y coordinate must be reversed before plortting but in this case are reversed directly in plotting    
        tmpSingleDisparities.append(tmpDisparities)
    

    # OverallDisparity is the computed disparity of the whole network in each
    # timewindow 
    # OverallDisparity[Timewindow]
    # will return a list containing the x and y disparities
    OverallDisparty.append((tmpDisparity/num_orientations))
    # SingleDisparities is the detected disparity of each orientation sensitive
    # Corolla in each timewindow 
    # SingleDisparities[Timewindow][Orientation][Corolla]
    # will return a list containing the x and y disparities
    SingleDisparities.append(tmpSingleDisparities)
    counter = counter + timeWindow
 
#%% 3D visualizer of the network activity
# Look at pyqtgraph documentation to learn how it's done.
# In the future using standard functions in the libraries (not stables at the moment of writing)
# it's kindly adviced

PlotName="Network decoded disparity 3D Scatterplot"
glViewWidget=gl.GLViewWidget()

## Add grid to the view
g = gl.GLGridItem()
g.setSize(x=stimuli_res_x,y=stimuli_res_y)
g.setSpacing(x=1,y=1)
g.setDepthValue(10)  # draw grid after surfaces since they may be translucent
g.translate(round(stimuli_res_x/2), round(stimuli_res_y/2), 0)
glViewWidget.addItem(g)

## Add axis to the view
a = gl.GLAxisItem()
a.setSize(x=stimuli_res_x/2,y=stimuli_res_y/2,z=200)
a.translate((stimuli_res_x-1)/2, (stimuli_res_y-1)/2, 0)
glViewWidget.addItem(a)


glViewWidget.show()
glViewWidget.setWindowTitle(PlotName)

#Dummy initialization
#X first and then Y this is how pyqtgraph expects points.
pos = np.zeros((stimuli_res_x, stimuli_res_y, 3))
pos[:, :, : 2] = np.mgrid[:stimuli_res_x, :stimuli_res_y].transpose(1, 2, 0)
pos = pos.reshape(stimuli_res_x*stimuli_res_y, 3)


#scattersplots

sp_fix = gl.GLScatterPlotItem(pos=pos, color=(1, 1, 1, .3), size=1,
                       pxMode=False)
spInL = gl.GLScatterPlotItem(pos=pos, color=(1, 1, 1, .3), size=1,
                       pxMode=False)
spInR = gl.GLScatterPlotItem(pos=pos, color=(1, 1, 1, .3), size=1,
                       pxMode=False)

spDisparity =  gl.GLScatterPlotItem(pos=pos, color=(1, 1, 1, .3), size=1,
                       pxMode=False)

spCorollasDisparities =  gl.GLScatterPlotItem(pos=pos, color=(1, 1, 1, .3), size=1,
                       pxMode=False)


glViewWidget.addItem(sp_fix)
glViewWidget.addItem(spInL)
glViewWidget.addItem(spInR)
glViewWidget.addItem(spDisparity)
glViewWidget.addItem(spCorollasDisparities)


glViewWidget.setCameraPosition(distance=200)
glViewWidget.pan(round(stimuli_res_x/2), round(stimuli_res_y/2), 0)

#adding some fixed point to the plot to clarify visualization (borders and center of the calculus window)
Perim = (window_size_x*2)+(window_size_y*2)-4 #this is the number of the cells in the perimeter
#the -4 prevent counting more than one cell per vertex
zoffs=0 #offset for the z value 

fixPos = np.zeros(shape=(Perim+1,3),dtype=np.float32)#+1 in order to add the center
fixPos[0,:] = [stimuli_res_x/2, stimuli_res_y/2, zoffs]
# Upper and Lower border
fixPos[1:window_size_x+1,:] = np.transpose([(np.arange(window_size_x))+int((stimuli_res_x-window_size_x+1)/2),
      np.ones(shape=(window_size_x),dtype=np.float32)*int((stimuli_res_y-window_size_y+1)/2),
      np.ones(shape=(window_size_x),dtype=np.float32)*zoffs])
startIndex = window_size_x+1
endIndex = (window_size_x*2)+1

fixPos[startIndex:endIndex, :] = np.transpose([(np.arange(window_size_x))+int((stimuli_res_x-window_size_x+1)/2),
      np.ones(shape=(window_size_x),dtype=np.float32)*int((stimuli_res_y+window_size_y-1)/2),
      np.ones(shape=(window_size_x),dtype=np.float32)*zoffs])
startIndex = endIndex
endIndex = startIndex+window_size_y-2
#Left and Right border
fixPos[startIndex:endIndex, :] = np.transpose([np.ones(shape=(window_size_y-2),
      dtype=np.float32)*int(fixPos[1,0]),
      (np.arange(1,window_size_y-1))+int((stimuli_res_y-window_size_y+1)/2),
      np.ones(shape=(window_size_y-2),dtype=np.float32)*zoffs])
startIndex = endIndex
endIndex = startIndex+window_size_y-2
fixPos[startIndex:endIndex, :] = np.transpose([np.ones(shape=(window_size_y-2),
      dtype=np.float32)*int(fixPos[window_size_x,0]),
      (np.arange(1,window_size_y-1))+int((stimuli_res_y-window_size_y+1)/2),
      np.ones(shape=(window_size_y-2),dtype=np.float32)*zoffs])

color = np.empty((len(fixPos), 4), dtype=np.float32)
color[:, 3] = 1
color[:, 0] = 1
color[:, 1] = 0
color[:, 2] = 0

fixPos[:,1] = stimuli_res_y -fixPos[:,1]
sp_fix.setData(pos=fixPos, color=color)



# Overall Disparity Vector
# obtained with the actvitiy of the whole population (global)
disparity_global_vector = []
NPoints=10
for k_timewindow in range(len(OverallDisparty)):
    tmpVector = np.zeros((NPoints,3))
    for k in range(NPoints):
        tmpVector[k,[0,1]] = ((OverallDisparty[k_timewindow]/5) + k*OverallDisparty[k_timewindow]/5) + np.asarray([stimuli_res_x/2 ,stimuli_res_y/2])
        tmpVector[k,2] = 0
    disparity_global_vector.append(tmpVector)

clock = 0
runtime = len(InposL)

# Signle Corolla's Disparity Vectors
# obtained with the actvitiy of each corolla (local)

disparity_local_vectors = []
NPoints = 10
for k_timewindow in range(len(OverallDisparty)):
    tmpsum = np.zeros((len(corollas_centers),2))
    for orient in range(num_orientations):
        tmpsum += SingleDisparities[k_timewindow][orient]
    tmpsum = tmpsum/num_orientations 
    allpoints = np.zeros((len(corollas_centers)*NPoints,3))
    for k in range(NPoints):
        allpoints[k*len(SingleDisparities[0][0]):k*len(SingleDisparities[0][0])+len(SingleDisparities[0][0]),
                  [0,1]] = ((tmpsum/10) + k*tmpsum/10) + corollas_centers + np.asarray([stimuli_res_x/2 ,stimuli_res_y/2]) - np.asarray([input_size_x/2 , input_size_y/2])
    allpoints[:,1] = stimuli_res_y - allpoints[:,1]
    allpoints[:,2] = 1    
    disparity_local_vectors.append(allpoints)


def update():
    global clock
    clock = clock*(clock<runtime) 
    print('/==============/')
    print(clock*timeWindow)
    # update surface positions and colors
    #### input population
    spInL.setData(pos=InposL[clock], color=IncolorsL[clock])
    spInR.setData(pos=InposR[clock], color=IncolorsR[clock])    
    ####Disparity decoded
    spDisparity.setData(pos=disparity_global_vector[clock], color=[0.,0.,1.,1.])
    spCorollasDisparities.setData(pos=disparity_local_vectors[clock], color=[0.,0.5,1.,1.])
    clock = clock + 1

t = QtCore.QTimer()
t.timeout.connect(update)
t.start(100)

#%% Stop visualisation from running

t.stop()


#%% Visualization - Corolla activity
# Plot that will show the activity of the Corollas centered in x0 and y0,
# the plot it is devided in each oriented corolla

# I am displaying Imem, therefore i need to normalize over the units,
# if you are not seeing anything, check a subset of disparity_neurons_corolla[orient].Imem[:]
# and set the right norm_cost
norm_const = mA

#I want to display each population with a different color
Popcolors = [[1, 0, 0], [89/255, 198/255, 118/255], [0, 0, 1], [247/255, 0, 255/255],
          [0, 1, 0], [1, 128/255, 0], [120/255, 120/255, 120/255], [0, 171/255, 1]]


pg.setConfigOptions(antialias=True)


labelStyle = {'color': '#FFF', 'font-size': '12pt'}
win = pg.GraphicsWindow(title='Activity Plot of Stereo Cell Corollas')
win.setWindowTitle('Mean activity of neurons centered in x='+str(x0)+' y='+str(y0)+' respective to L plane')  
plots = []
scatters = []
Trigger=0
Relative_max = []
for k in range(num_orientations):  
    # trigger to switch to the next row in the plot
    if(k+1 > num_orientations/2) and not Trigger:
        win.nextRow()
        Trigger=1
    plots.append(win.addPlot(title=str(k)+"/8"+"π Sensitive Corolla"))
    scatters.append(pg.ScatterPlotItem(pxMode=False, size=1, brush=pg.mkBrush(Popcolors[k])))
    plots[k].addItem(scatters[k])
    pos = neuron_centers[ind,2:4]
    scatters[k].setData(pos=pos)
    scatters[k].getViewBox().invertY(True)
    Relative_max.append(np.max(disparity_neurons_corolla[k].Imem))

up_start=0
up_end=len(disparity_neurons_corolla[0].t)    


up_index=up_start
Relative_max[4]
spots = []
def update():
    global up_index, spots
    up_index += 1
    if(up_index>up_end):
        up_index = up_start-1
    for k in range(num_orientations):
        # z is the activity of the subplot
        z = disparity_neurons_corolla[k].Imem[:,up_index]/(np.max(Relative_max)*norm_const)
        spots.clear()
        spots = [{'pos': (pos[n][0], pos[n][1]), 'brush': pg.mkBrush((int(Popcolors[k][0]*z[n]), int(Popcolors[k][1]*z[n]), int(Popcolors[k][2]*z[n])))} for n in range(len(ind))]
        #adding the artificial disparity, comment if you want to visualize recorded data
        spots.append({'pos': (xoff+x0,yoff+y0), 'brush': pg.mkBrush(255, 255, 255,), 'symbol': '+'})         
        scatters[k].setData(spots)
    print(up_index)

t = QtCore.QTimer()
t.timeout.connect(update)
t.start(1)    

#%% Stop visualisation from running

t.stop()


#%% Visualization - Single Corolla of oriented disparity and inhibitory neurons, membrane Current
# This plot serves the scope of assest visually the ability of inhibitory neurons
# to reduce false matching.
# If the inhibitory neurons are firing but they are not stopping disparity neurons
# something is wrong.

#I select the Oriented layer that i want to display and i encode distance from the center with colors
Orient=4

#I want to display each population with a different color
Popcolors = [(255, 0, 0), (89, 198, 118), (0, 0, 255), (247, 0, 255),
          (0, 255, 0), (255, 128, 0), (120, 120, 120), (0, 171, 255)]

labelStyle = {'color': '#FFF', 'font-size': '12pt'}
win2 = pg.GraphicsWindow(title='Stereo and Inhibitory cells state')
win2.setWindowTitle('Stereo and Inhibitory cells state')

p1 = win2.addPlot(title='Stereo cells')
win2.nextRow()
p2 = win2.addPlot(title='Inhibitory cells L')
win2.nextRow()
p3 = win2.addPlot(title='Inhibitory cells R')



colorlevels_x=int(255/(np.max(np.dot(np.abs(pos[:,0]-x0),1/(spacing+1)))+1))
colorlevels_y=int(255/(np.max(np.dot(np.abs(pos[:,1]-y0),1/(spacing+1)))+1))
#I encode relative distance with 255 int values, it won't work with more than 500 "petals"
Distance_x=np.zeros(len(ind))
Distance_y=np.zeros(len(ind))
for k in range(len(ind)):
    Distance_x[k] = 255- int(np.dot(np.abs(pos[k,0]-x0),1/(spacing+1)))*colorlevels_x
    Distance_y[k] = 255- int(np.dot(np.abs(pos[k,1]-y0),1/(spacing+1)))*colorlevels_y
    
Beginshw=100
Endshw=2000

for k, data in enumerate(np.asarray(disparity_neurons_corolla[Orient].Imem)):
    color = (Distance_x[k],Distance_y[k],255)
    p1.plot(x=np.asarray(disparity_neurons_corolla[Orient].t[Beginshw:Endshw] / ms), y=data[Beginshw:Endshw],
    pen=pg.mkPen(color, width=2))
    
for k, data in enumerate(np.asarray(inhibitory_neurons_L_corolla[Orient].Imem)):
    color = (Distance_x[k],Distance_y[k],255)
    p2.plot(x=np.asarray(inhibitory_neurons_L_corolla[Orient].t[Beginshw:Endshw] / ms), y=data[Beginshw:Endshw],
    pen=pg.mkPen(color, width=2))
    
for k, data in enumerate(np.asarray(inhibitory_neurons_R_corolla[Orient].Imem)):
    color = (Distance_x[k],Distance_y[k],255)
    p3.plot(x=np.asarray(inhibitory_neurons_R_corolla[Orient].t[Beginshw:Endshw] / ms), y=data[Beginshw:Endshw],
    pen=pg.mkPen(color, width=2))  

# Normalizing and scaling the plots    
yMax= np.max([np.max(disparity_neurons_corolla[Orient].Imem),np.max(inhibitory_neurons_L_corolla[Orient].Imem),np.max(inhibitory_neurons_R_corolla[Orient].Imem)])
p1.setXRange(up_start, up_end, padding=0)
p2.setXRange(up_start, up_end, padding=0)
p3.setXRange(up_start, up_end, padding=0)
p1.setYRange(0, yMax, padding=0)
p2.setYRange(0, yMax, padding=0)
p3.setYRange(0, yMax, padding=0)
p1.setLabel('left', "Imem", units='A', **labelStyle)
p1.setLabel('bottom', "Time (ms)", **labelStyle)
p2.setLabel('left', "Imem", units='A', **labelStyle)
p2.setLabel('bottom', "Time (ms)", **labelStyle)
p3.setLabel('left', "Imem", units='A', **labelStyle)
p3.setLabel('bottom', "Time (ms)", **labelStyle)

#%% Firing analysis
# This script computes the 
zero_disparity_pos =  ((neuron_centers[:,0]==neuron_centers[:,2])*1+(neuron_centers[:,1]==neuron_centers[:,3])*1)==2
zero_disparity_ind = np.arange(num_stereo_cells)[zero_disparity_pos]
Correct_matches = 0
for k in range(len(stereo_spikes[4].i)):
    for kk in range(len(zero_disparity_ind)):
        if(stereo_spikes[4].i[k]==zero_disparity_ind[kk]):
            Correct_matches += 1
False_matches = len(stereo_spikes[4].i) - Correct_matches

Succes_rate = Correct_matches/False_matches *100
print(Succes_rate)
print(len(stereo_spikes[4].i))


#%% Save disparity_neurons_activity
# This is basically a script to export the result in .mat

OutMat = {}  #An output specifically designed for Matlab
for k in range(num_orientations):
    OutMat['OnPop'+str(k)] =  disparity_neurons_activity[str(k)]
    OutMat['OnPop'+str(k)+'t'] =  disparity_neurons_activity[str(k)+'t']
#    OutMat['OffPop'+str(k)] =  disparity_neurons_activity[str(k)]
#    OutMat['OffPop'+str(k)+'t'] =  disparity_neurons_activity[str(k)+'t']

OutMat['start'] = disparity_neurons_activity['start']
OutMat['end'] = disparity_neurons_activity['end']




sio.savemat(parent_folder+'/Data/Parvo_Matlab_Out/Input_On_Moving_Bar.mat', OutMat)



