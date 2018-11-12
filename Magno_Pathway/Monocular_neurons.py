#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 14:52:49 2017

@author: marcorax93

This script it has been developed on spyder, it is devided in sometimes
indipendent subsections, marked with #%%, thus do not run the entire script as whole

Using spyder is advised, also remember that the plots make use of pyqtgraph, 
for its speed on intensive 3d plotting, for this reason remember to disable
matplotlib default support on Spyder (Tools/preferences/IPython console/Graphics)
 
"""

# General computation and plotting
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import numpy as np
import pickle
import time
from pathlib import Path
from brian2 import ms, mV, pA, nS, nA, pF, us, volt, second, prefs,\
    SpikeMonitor, StateMonitor, figure, plot, show, xlabel, ylabel,\
    seed, xlim, ylim, subplot, network_operation, TimedArray,\
    defaultclock, SpikeGeneratorGroup, asarray, pamp, set_device, device
# Teili lib
from teili.core.groups import Neurons, Connections
from teili import teiliNetwork
from teili.models.neuron_models import DPI 
from teili.models.synapse_models import DPISyn
from teili.tools.converter import aedat2numpy, dvs2ind
# My libs
from Tools.Stimuli import rotating_bar
from Tools.Davisloading import  eventsToPy, fromNPtoBrian
from Tools.Connectivity import connect_spike_gen, gabor_spike_gen



# Settings for c++ code building, it helps to speed up the computation a lot in 
# long simulations 
prefs.codegen.target = "numpy"
set_device('cpp_standalone', directory='C++ Buids/Mono_R_Off_Moving_Bar')
prefs.devices.cpp_standalone.extra_make_args_unix = ['-j 16'] #The total number of file being built in parallel at once (this value depend on the number of cores, schedulers and so on) tweak it to find the best value on your machine
prefs.devices.cpp_standalone.openmp_threads = 8 #The maximum number of threads the run will feed 

#The default clock of the simulation
defaultclock.dt = 1 * ms

#Important to find the Path to the Data folder
parent_folder=str(Path().resolve().parent)

#%% Net Parameters - Recorded Stimuli
#run this only if you want to load and run recorded stimuli

num_orientations = 8 # number of subpopulations of neurons sensitive to orientations
theta=np.linspace(0., (7/8)*np.pi, num=num_orientations) # Orientation angle array,
                                                         # defining the angle which each subpopulation
                                                         # will be sensitive to 
kernel_size=31 # linear size of the squared connectivity kernel used to connect
               # simple cells to the image plane aka receptive field size
               # it has to be an odd number

kernel_borders = (kernel_size-1)/2

# The kernel is gabor like and the needed parameters are these 
sigma_x=8.4  
sigma_y=7.
freq=11.

scal_weight=1# this is a scaling weight applied to all the synapses of a connectivity kernel 

input_size_x=240  # lateral resolutions of the image plane
input_size_y=180


net_size_x = 40 # lateral dimensions of the network, it can be only smaller than the image plane
net_size_y = 80 # before setting these values check if you'll end up having neurons
                # on the borders of the image (in other words if their connectivity kernels)
                # won't connect to non existent neurons outside the image plane.
                # check if 
                # net_size_x < input_size_x - kernel_borders
                # net_size_y < input_size_y - kernel_borders

offx = int((input_size_x-net_size_x)/2  + 0 ) # basically where do you want to center the network in respect to image coordinates
offy = int((input_size_y-net_size_y)/2  + 0 )  # don't get out of borders with the connectivity kernels! 

image_res = input_size_x*input_size_y

num_neurons_sub_pop = net_size_x*net_size_y # the number of neurons in each sub population

total_num_neurons = num_neurons_sub_pop*num_orientations


#%% Net Parameters - Artifiial Stimuli
#run this only if you want to load and run artificial stimuli

num_orientations = 8 # number of subpopulations of neurons sensitive to orientations
theta=np.linspace(0., (7/8)*np.pi, num=num_orientations) # Orientation angle array,
                                                         # defining the angle which each subpopulation
                                                         # will be sensitive to 
kernel_size=11 # linear size of the squared connectivity kernel used to connect
               # simple cells to the image plane aka receptive field size
               # it has to be an odd number

kernel_borders = (kernel_size-1)/2

# The kernel is gabor like and the needed parameters are these 
sigma_x=3. 
sigma_y=2.5
freq=4.

scal_weight=60# this is a scaling weight applied to all the synapses of a connectivity kernel 

input_size_x=240  # lateral resolutions of the image plane
input_size_y=180


net_size_x = 40 # lateral dimensions of the network, it can be only smaller than the image plane
net_size_y = 40 # before setting these values check if you'll end up having neurons
                # on the borders of the image (in other words if their connectivity kernels)
                # won't connect to non existent neurons outside the image plane.
                # check if 
                # net_size_x < input_size_x - kernel_borders
                # net_size_y < input_size_y - kernel_borders

offx = int((input_size_x-net_size_x)/2  + 0 ) # basically where do you want to center the network in respect to image coordinates
offy = int((input_size_y-net_size_y)/2  + 0 )  # don't get out of borders with the connectivity kernels! 

image_res = input_size_x*input_size_y

num_neurons_sub_pop = net_size_x*net_size_y # the number of neurons in each sub population

total_num_neurons = num_neurons_sub_pop*num_orientations

 
#%% DAVIS data extraction and saving - aedat 3.1, 3.0 - Recorded stimuli
# Run this part only once to save everything in PopulationData
numpyevents=eventsToPy(filename=parent_folder+'/Data/DAVIS_Recordings/Events_R_Moving_Bar-2018_03_06_17_04_05.aedat', xdim=input_size_x, ydim=input_size_y)


#final conversion and data saving
indInp2 = fromNPtoBrian(Events=numpyevents, xRes=input_size_x, uniqueflag=1, milliflag=1)


# Saving for later (useful if you want to run simulation on the same data multiple time )
outputFile = parent_folder+'/Data/Magno_Population_Data/R_Moving_Bar_input.data'
fw = open(outputFile, 'wb')
pickle.dump(indInp2, fw)
fw.close()





#%% DAVIS data loading, chopping and variable renaming - Recorded stimuli

# Loading data 
inputFile = parent_folder+'/Data/Magno_Population_Data/R_Moving_Bar_input.data' 
fd = open(inputFile, 'rb')
indInp2 = pickle.load(fd)

# Selecting begin and of the dataset, useful for focusing on a particular time window
start = 0*second
end = np.max(indInp2[0][1])*ms

# Chopping data and renaming variables(explicit separation of ON and OFF events)
indOn=indInp2[0][0][(indInp2[0][1][:]*us >= start)*(indInp2[0][1][:]*us <= end)]
indOff=indInp2[1][0][(indInp2[1][1][:]*us >= start)*(indInp2[1][1][:]*us <= end)]
timeOn=indInp2[0][1][(indInp2[0][1][:]*us >= start)*(indInp2[0][1][:]*us <= end)]
timeOff=indInp2[1][1][(indInp2[1][1][:]*us >= start)*(indInp2[1][1][:]*us <= end)]
timeOn=(timeOn*ms)-start
timeOff=(timeOff*ms)-start


# select if you want to process ON or OFF events 

retina = SpikeGeneratorGroup(image_res, indices=indOn, times=timeOn, name='gPre')

#retina = SpikeGeneratorGroup(image_res, indices=indOff, times=timeOff, name='gPre')



#%% Artificial stimuli 
# Rotating bar with known angles, used to assest the net ability to encode orientation

length = 15
center_offx = 0
center_offy = 0
angle_step = 40
ts_offset = 10 * defaultclock.dt

retina, events = rotating_bar(length,  input_size_x, input_size_y, center_offx,
                 center_offy, angle_step, ts_offset)

angles_bar_degrees = np.linspace(0, 360, angle_step)

start = 0 * ms  
end = events[2][-1:]
end = end[0]*ts_offset


#%% Network creation

start_time = time.time()
Net = teiliNetwork()

# Model (build the equation string with the equation builder from tili)
# the number of inputs need to be equal to the number of populations 
# this is done to ensure that each population will have different inputs names
neuron_model = DPI(num_inputs=num_orientations)

# x and y cooridinates for each monocular neuron in the net, in other words,
# the retinal coordinates where each monocular connectivity kernel will be 
# centered in 
spatial_positions_x = (np.arange(num_neurons_sub_pop) % net_size_x) + offx
spatial_positions_y =  (np.arange(num_neurons_sub_pop) // net_size_x) + offy

# Neurons
populations = []

#(I am using the same model for all the neurons)
for k in range(num_orientations):       
    populations.append(Neurons(N=num_neurons_sub_pop, name="simple_neurons_pop_"+str(k), equation_builder=neuron_model))
    # Set spatial position for each neuron in the subpopublation
    populations[k].x = spatial_positions_x
    populations[k].y = spatial_positions_y

# Model (build the equation string with the equation builder from tili)
synapse_model = DPISyn()

# The synapses that will connect each subpopulation of simple cells to the "retinal or image plane"
synapses = []

# At the moment the only way to load functions with addictional parameters
# it's through string sobstitution, I don't know if it's a tili or brian2 problem
synapses_params = {'input_size_x': input_size_x,
                   'input_size_y': input_size_y,
                   'kernel_size': kernel_size,
                   'scal_weight': scal_weight,
                   'sigma_x': sigma_x,
                   'sigma_y': sigma_y,
                   'freq': freq}

for k in range(num_orientations):       
    synapses.append(Connections(retina, populations[k],
                     name="InpSyn_pop_"+str(k), equation_builder=synapse_model))
    # x and y are casted as integer because the implemented connectivity kernels are supposed to work with indeces, not with positions
    synapses[k].connect('connect_spike_gen(i, int(x_post), int(y_post), {input_size_x}, {input_size_y}, {kernel_size})'.format(**synapses_params))
    synapses_params['theta'] = theta[k]
    # gabor retina is not added to the namespace at the time of the creation of neuron and synapses(weird)
    synapses[k].namespace['gabor_spike_gen'] = gabor_spike_gen
    synapses[k].weight = '{scal_weight}*gabor_spike_gen(i, int(x_post), int(y_post), {input_size_x}, {input_size_y}, {theta}, {sigma_x}, {sigma_y}, {freq})'.format(**synapses_params)

#Monitors
    
# Spiking Activity

# Monitor the spiking activity of the input  
spikemon_retina = SpikeMonitor(retina, name='spikemon_retina')    

# Monitor the spiking activity of the populations of simple neurons  
neurons_spikes = []
for k in range(num_orientations):       
    neurons_spikes.append(SpikeMonitor(populations[k], name='spikemon_pop_'+str(k)))

# Monitor current activity of a single neuron per population
# Useful to debug the model or some parameters
ind=num_neurons_sub_pop/2
neurons_states = []
for k in range(num_orientations):       
    neurons_states.append(StateMonitor(populations[k], variables=["Iin", "Imem"], record=ind, name='statemon_pop_'+str(k)))
    
# Monitor current activity of a single synapse per population
# Useful to debug the model or some parameters
ind=(num_neurons_sub_pop*kernel_size**2)/2
synapse_states = []
for k in range(num_orientations):       
    synapse_states.append(StateMonitor(synapses[k], variables=["Ie_syn", "Ii_syn"], record=ind, name='statemon_syn_'+str(k)))

Net.add(retina, populations, 
        synapses, spikemon_retina, neurons_spikes, neurons_states, synapse_states)

# Tili gives the possibilty to change some parameters called standalone parameters
# without rebuilding the network, the support with subgroups at the moment of writing
# it's sadly still experimental

#Net.add_standalone_params(scal_weight=scal_weight)

Net.run(end-start)
print("Network computation time, I hope the coffe was good at least: " + str(time.time() - start_time))



#%% Connectivety plot 3D

colors = [(255, 0, 0), (89, 198, 118), (0, 0, 255), (247, 0, 255),
          (0, 255, 0), (255, 128, 0), (120, 120, 120), (0, 171, 255)]

# select the orientation you want to plot
orient = 0
gabor_orientation = theta[orient]

## Create a GL View widget to display data
w = gl.GLViewWidget()
w.setWindowTitle('3D Gabor kernel theta:'+str(gabor_orientation/np.pi)+'Pi')

# set pan 
w.pan((kernel_size-1)/2,(kernel_size-1)/2,0)

scaling_fact = 2
z = scaling_fact*gabor_spike_gen(range(kernel_size**2), (kernel_size-1)/2, (kernel_size-1)/2,  kernel_size, kernel_size, gabor_orientation, sigma_x, sigma_y, freq).reshape(kernel_size,kernel_size)
# NB I flipped upsidedown the matrix z because the reference for plotting is different from the population's
p1 = gl.GLSurfacePlotItem(z=np.flipud(np.transpose(z)), shader='heightColor', color=(0.5, 0.5, 1, 1))
w.addItem(p1)


## Add a grid to the view
g = gl.GLGridItem()
g.setSize(x=kernel_size,y=kernel_size)
g.setSpacing(x=1,y=1)
g.setDepthValue(10)  # draw grid after surfaces since they may be translucent
g.translate((kernel_size-1)/2, (kernel_size-1)/2, 0)
w.addItem(g)
## Add axis to the view
a = gl.GLAxisItem()
a.setSize(x=kernel_size/2,y=kernel_size/2,z=np.max(z)*2)
a.translate((kernel_size-1)/2, (kernel_size-1)/2, 0)
w.addItem(a)
w.show()

#%% Activity

colors = [(255, 0, 0), (89, 198, 118), (0, 0, 255), (247, 0, 255),
          (0, 255, 0), (255, 128, 0), (120, 120, 120), (0, 171, 255)]

#Plotting this graph can kill the kernel if the actvity is high and the recorded
#neuron and synapses are too many, to preventi this from happening , plot only 
# a small window of activity

#Please remember that the begin and end here plotted are the number of the steps
# of the default clock at the beginning of this file
begin_2D_plot = 0 * defaultclock.dt
end_2D_plot = 500 * defaultclock.dt

labelStyle = {'color': '#FFF', 'font-size': '12pt'}
win2 = pg.GraphicsWindow(title='Simple cells response M Pathway')
win2.setWindowTitle('Simple cells response')

p1 = win2.addPlot(title="Retina activity")
win2.nextRow()
p2 = win2.addPlot(title="Simple neurons activity")
win2.nextRow()
p3 = win2.addPlot(title="Input synapses states")
win2.nextRow()
p4 = win2.addPlot(title='Simple neuron states')



## Retina spiking activity
timetmp = (spikemon_retina.t>=begin_2D_plot)*(spikemon_retina.t<=end_2D_plot) # the indeces of all events in the selected time window
p1.plot(x=np.asarray(spikemon_retina.t[timetmp] / ms), y=np.asarray(spikemon_retina.i[timetmp]),
        pen=None, symbol='o', symbolPen=None,
        symbolSize=7, symbolBrush=(255, 255, 255))

## Simple neuron spiking activity
for k in range(num_orientations):
    timetmp = (neurons_spikes[k].t>=begin_2D_plot)*(neurons_spikes[k].t<=end_2D_plot) # the indeces of all events in the selected time window
    p2.plot(x=np.asarray(neurons_spikes[k].t[timetmp] / ms), y=np.asarray(neurons_spikes[k].i[timetmp]),
        pen=None, symbol='o', symbolPen=None,
        symbolSize=7, symbolBrush=colors[k])
    
    

## Synapses states
p3.addLegend()

# Excitatory
for k in range(num_orientations):
    nI = np.shape(synapse_states[k].Ie_syn)[0]
    for j, data in enumerate(np.asarray(synapse_states[k].Ie_syn)):
        timetmp = (synapse_states[k].t>=begin_2D_plot)*(synapse_states[k].t<=end_2D_plot) # the indeces of all events in the selected time window
        name = 'excitatory_synapse_orientation_{}'.format(k)
        color=  (colors[1][0], colors[2][1], j*int(np.floor(255/np.shape(synapse_states[k].Ii_syn)[0])))
        print (name, "ind"+str(j), color)
        p3.plot(x=np.asarray(synapse_states[k].t[timetmp] / ms), y=data[timetmp],
                pen=pg.mkPen(color, width=1), name=name)


# Inibitory
for k in range(num_orientations):
    nI = np.shape(synapse_states[k].Ii_syn)[0]
    for j, data in enumerate(np.asarray(synapse_states[k].Ii_syn)):
        timetmp = (synapse_states[k].t>=begin_2D_plot)*(synapse_states[k].t<=end_2D_plot) # the indeces of all events in the selected time window
        name = 'inhibitory_synapse_orientation_{}'.format(k)
        color=  (colors[3][0], colors[3][1], j*int(np.floor(255/np.shape(synapse_states[k].Ii_syn)[0])))
        print (name, "ind"+str(j), color)
        p3.plot(x=np.asarray(synapse_states[k].t[timetmp] / ms), y=-data[timetmp],
                pen=pg.mkPen(color, width=1), name=name)




## Simple neurons states
for k in range(num_orientations):
    for j, data in enumerate(np.asarray(neurons_states[k].Imem)):
        timetmp = (neurons_states[k].t>=begin_2D_plot)*(neurons_states[k].t<=end_2D_plot) # the indeces of all events in the selected time window
        p4.plot(x=np.asarray(neurons_states[k].t[timetmp] / ms), y=data[timetmp],
        pen=pg.mkPen(colors[k], width=2))



p1.setLabel('left', "Neuron ID", **labelStyle)
p1.setLabel('bottom', "Time (ms)", **labelStyle)
p2.setLabel('left', "Neuron ID", **labelStyle)
p2.setLabel('bottom', "Time (ms)", **labelStyle)
p3.setLabel('left', "Synaptic current Isyn", units='A', **labelStyle)
p3.setLabel('bottom', "Time (ms)", **labelStyle)
p4.setLabel('left', "Membrane current Imem", units="A", **labelStyle)
p4.setLabel('bottom', "Time (ms)", **labelStyle)

minx=p4.viewRange()[0][0]
maxx=p4.viewRange()[0][1]

b = QtGui.QFont("Sans Serif", 10)
p1.getAxis('bottom').tickFont = b
p1.getAxis('left').tickFont = b
p1.setXRange(minx, maxx, padding=0)
p2.getAxis('bottom').tickFont = b
p2.getAxis('left').tickFont = b
p2.setXRange(minx, maxx, padding=0)
p3.getAxis('bottom').tickFont = b
p3.getAxis('left').tickFont = b
p3.setXRange(minx, maxx, padding=0)
p4.getAxis('bottom').tickFont = b
p4.getAxis('left').tickFont = b





#%% Frametizer 
# In order to display an animated output of the network activity
# I need to produce frames of events.
# It also computes a vector displaying the direction of the stimuli decoded 
# in the window. NB it will take a while! Have a second coffe


timeWindow = 10 * defaultclock.dt # The width in time used to build a frame out of single events

counter = 0 * defaultclock.dt

indexarrayIn=spikemon_retina.i[:]
timearrayIn=spikemon_retina.t[:]

#I want to display each population with a different color
Popcolors = [[1, 0, 0], [89/255, 198/255, 118/255], [0, 0, 1], [247/255, 0, 255/255],
          [0, 1, 0], [1, 128/255, 0], [120/255, 120/255, 120/255], [0, 171/255, 1]]

Inpos = []
Outpos = []
Incolors = []
Outcolors = []
DirectionStimuli = [] # the resulting direction of the stimuli of the whole net
while ((end-start)>counter):
    timetmp = (timearrayIn>=counter)*(timearrayIn<counter+timeWindow) # the indeces of all events in this current time window
    z = np.zeros(len(timearrayIn[timetmp==1])) # z coordinate
    tmpindex=np.unravel_index(indexarrayIn[timetmp==1], (input_size_y, input_size_x))    
    # the x coordinates are already right, the y coordinate have to be mirrored. The coordinate references are different from the frame to the plot 
    Inpos.append(np.transpose(np.asarray([tmpindex[:][1], -tmpindex[:][0] + input_size_y, z], dtype=np.float32 )))
    
    color = np.empty((len(z), 4), dtype=np.float32)
    color[:, 3] = 0.9
    color[:, 0] = np.clip(5.0, 0, 1)
    color[:, 1] = np.clip(3, 0, 1)
    color[:, 2] = np.clip(1.0, 0, 1)
    Incolors.append(color)
    tmpOutpos = [] 
    tmpOutcolors = []
    tmpPopulationVectors = [] 
    tmpXcoord = 0
    tmpYcoord = 0
    for k in range(num_orientations):
        timetmp = (neurons_spikes[k].t>=counter)*(neurons_spikes[k].t<counter+timeWindow)
        z = abs(np.sin(((neurons_spikes[k].t[timetmp==1]-counter)/us)+np.pi))+(k*4)
        tmpindex=np.unravel_index(neurons_spikes[k].i[timetmp==1], (net_size_y, net_size_x))
        tmpOutpos.append(np.transpose(np.asarray([tmpindex[:][1] + offx,
                                   (-tmpindex[:][0] + net_size_y + offy),
                                   z+2], dtype=np.float32 )))
        color = np.empty((len(z), 4), dtype=np.float32)
        color[:, 3] = 0.9
        color[:, 0] = np.clip(z * Popcolors[k][0], 0, 1)
        color[:, 1] = np.clip(z * Popcolors[k][1], 0, 1)
        color[:, 2] = np.clip(z * Popcolors[k][2], 0, 1)
        tmpOutcolors.append(color)
        tmpXcoord = tmpXcoord + len(tmpOutpos[k]) * np.cos((np.pi*k/8)) 
        tmpYcoord = tmpYcoord + len(tmpOutpos[k]) * np.sin((np.pi*k/8))
    #tmpXcoord = tmpXcoord + len(tmpOutpos[0]) * np.cos((np.pi*0/8)) 
    #tmpYcoord = tmpYcoord + len(tmpOutpos[0]) * np.sin((np.pi*0/8))
    tmpPopulationVectors.append([tmpXcoord, tmpYcoord]) 
    Outpos.append(tmpOutpos)   
    Outcolors.append(tmpOutcolors)
    
    if (sum(tmpPopulationVectors[0]) != 0) :
        DirectionStimuli.append(np.angle(tmpXcoord+tmpYcoord*1j, deg=True))
    else:   
        DirectionStimuli.append(0)
    counter = counter + timeWindow

#%% Angle decoder-artificial stimuli 
#Small plot that shows the angle of the artificial stimuli, and the angle
# obtained by population decoding in the network, can be used to check if the 
# hypercolumn organization is working.
    
win = pg.GraphicsWindow(title="Window decoded angle ")
win.resize(1000,600)
win.setWindowTitle("Window decoded angle")
labelStyle = {'color': '#FFF', 'font-size': '12pt'}

# Enable antialiasing for prettier plots
pg.setConfigOptions(antialias=True)

p1 = win.addPlot(title="Decoded Angle", y=DirectionStimuli, pen=(200,200,200), symbolBrush=(255,0,0), symbolPen='w')
p1.plot(y=(angles_bar_degrees), pen=(200,100,200), symbolBrush=(0,255,0), symbolPen='w')
p1.setLabel('left', "Angle in degrees", **labelStyle)   
p1.setLabel('bottom', "Time bin", **labelStyle)   


    
#%% 3D visualizer of the network activity
# Look at pyqtgraph documentation to learn how it's done.
# In the future using standard functions in the libraries (not stables at the moment of writing)
# it's kindly adviced
PlotName="Network spiking actvity 3D Scatterplot"
glViewWidget=gl.GLViewWidget()

# Add grid to the view
g = gl.GLGridItem()
g.setSize(x=input_size_x, y=input_size_y)
g.setSpacing(x=1,y=1)
g.setDepthValue(10)  # draw grid after surfaces since they may be translucent
g.translate(round(input_size_x/2), round(input_size_y/2), 0)
glViewWidget.addItem(g)

# Add axis to the view
a = gl.GLAxisItem()
a.setSize(x=input_size_x/2,y=input_size_y/2,z=200)
a.translate((input_size_x-1)/2, (input_size_y-1)/2, 0)
glViewWidget.addItem(a)


glViewWidget.show()
glViewWidget.setWindowTitle(PlotName)

# Dummy initialization
# Note that pyqtgraph expects x first and then y coordinates
pos = np.zeros((input_size_x, input_size_y, 3))
pos[:, :, : 2] = np.mgrid[:input_size_x, :input_size_y].transpose(1, 2, 0)
pos = pos.reshape(input_size_x*input_size_y, 3)


# scattersplots 
sp_fix = gl.GLScatterPlotItem(pos=pos, color=(1, 1, 1, .3), size=1,
                       pxMode=False)
spIn = gl.GLScatterPlotItem(pos=pos, color=(1, 1, 1, .3), size=1,
                       pxMode=False)
glViewWidget.addItem(sp_fix)
glViewWidget.addItem(spIn)

simpleScatterPlot = []
for k in range(num_orientations):    
    simpleScatterPlot.append(gl.GLScatterPlotItem(pos=pos, color=(1, 1, 1, .3), size=1,
                       pxMode=False))
    glViewWidget.addItem(simpleScatterPlot[k])


glViewWidget.setCameraPosition(distance=200)
glViewWidget.pan(round(input_size_x/2), round(input_size_y/2), 0)

startv= start  
runtime = end-start

clock = startv

#adding some fixed point to the plot to clarify visualization (borders and center of the calculus window)
Perim = (net_size_x*2)+(net_size_y*2)-4 #this is the number of the cells in the perimeter
# the -4 prevent counting more than one cell per vertex
zoffs=0 #offset for the z value 

fixPos = np.zeros(shape=(Perim+1,3),dtype=np.float32)#+1 in order to add the center
fixPos[0,:] = [input_size_x/2, input_size_y/2, zoffs]
# Upper and Lower border
fixPos[1:net_size_x+1,:] = np.transpose([(np.arange(net_size_x))+offx,
      np.ones(shape=(net_size_x),dtype=np.float32)*offy,
      np.ones(shape=(net_size_x),dtype=np.float32)*zoffs])
startIndex = net_size_x+1
endIndex = (net_size_x*2)+1

fixPos[startIndex:endIndex, :] = np.transpose([(np.arange(net_size_x))+offx,
      np.ones(shape=(net_size_x),dtype=np.float32)*(offy+net_size_y-2),
      np.ones(shape=(net_size_x),dtype=np.float32)*zoffs])
startIndex = endIndex
endIndex = startIndex+net_size_y-2
# Left and Right borders
fixPos[startIndex:endIndex, :] = np.transpose([np.ones(shape=(net_size_y-2),
      dtype=np.float32)*int(fixPos[1,0]),
      (np.arange(1,net_size_y-1))+offy,
      np.ones(shape=(net_size_y-2),dtype=np.float32)*zoffs])
startIndex = endIndex
endIndex = startIndex+net_size_y-2
fixPos[startIndex:endIndex, :] = np.transpose([np.ones(shape=(net_size_y-2),
      dtype=np.float32)*int(fixPos[net_size_x,0]),
      (np.arange(1,net_size_y-1))+offy,
      np.ones(shape=(net_size_y-2),dtype=np.float32)*zoffs])

color = np.empty((len(fixPos), 4), dtype=np.float32)
color[:, 3] = 1
color[:, 0] = 1
color[:, 1] = 0
color[:, 2] = 0
sp_fix.setData(pos=fixPos, color=color)


clock = 0
runtime = len(Inpos)
# This is the function that will run in background
def update():
    global clock
    clock = clock*(clock<runtime) 
    print('/==============/')
    print(clock*timeWindow)
    # update surface positions and colors
    #### input population
    spIn.setData(pos=Inpos[clock], color=Incolors[clock])
    ####simple cells populations
    for k in range(num_orientations):
        simpleScatterPlot[k].setData(pos=Outpos[clock][k], color=Outcolors[clock][k])
    clock = clock + 1

t = QtCore.QTimer()
t.timeout.connect(update)
t.start(100)

#%% Stop visualisation from running

t.stop()


#%% Save network spiking activity

pop_data = {}
for k in range(num_orientations):    
    pop_data[str(k)] = neurons_spikes[k].i[:]
    pop_data[str(k)+'t'] = neurons_spikes[k].t[:]

pop_data['start'] = start
pop_data['end'] = end

outputFile = parent_folder+'/Data/Magno_Population_Data/Rotating_Artificial_Bar.data'
fw = open(outputFile, 'wb')
pickle.dump(pop_data, fw)
fw.close()


