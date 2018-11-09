#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 6 14:52:49 2018

@author: marcorax93

This file contain the artificial stimuli that can be used to test the network
 
"""

import numpy as np
import operator
from brian2 import SpikeGeneratorGroup,ms

def dda_round(x):
        return (x + 0.5).astype(int)

def xy2ind(x, y, xdim):
    """Given a pair of x, y (pixel) coordinates this function
    will return an index that correspond to a flattened pixel array

    Args:
        x (int, required): x-coordinate
        y (int, required): y-coordinate
        ydim (int): Row length

    Returns:
        ind (int): Converted index (e.g. flattened array)
    """
    if isinstance(x, np.ndarray):
        return x + (y * xdim)
    else:
        return int(x) + int(y) * xdim
    
def rotating_bar(length=10,  input_size_x=240, input_size_y=180, offx=0,
                 offy=0, angle_step=10, ts_offset=10*ms):
    """
    This function returns a single spikegenerator group (Brian object),
    and the events related to it.
    It produces a stimuli of a rotating bar, ideal for testing orientation 
    selective populations
    
    Args:
        length (int): `length` of the bar in pixel.
        input_size_x(int) : size of the input population x coordinate
                            by default set as the DAVIS resolution
        input_size_y(int) : size of the input population y coordinate
                            by default set as the DAVIS resolution
        offx(int): horizontal offset of the stimuli center (0, or centered in the 
                                                        middle by default)
        offy(int): vertical offset of the stimuli center (0, or centered in the 
                                                        middle by default)
        angle step(int) : the number of steps that the bar will use in order
                            to rotate starting from 0 to 2*pi
        ts_offset (seconds_brian): time between each consecutive fire of the population
                                        (time occurring between each bar position)
                                        It's only applied in the spike generator,
                                        the event timestamps don't have units!
    Returns:
        SpikeGenerator obj: Brian2 objects which holds the spiketimes as well
            as the respective neuron indices
        events list: A list containing the events as x,y coordinates, polatirty
                    and timestamps
    """
    x_coord = []
    y_coord = []
    pol = []
    ts = []
    center = (input_size_x/2 + offx, input_size_y/2 + offy)
    angles = np.linspace(0, np.pi*2, angle_step)
    for i, cAngle in enumerate(angles):
        endy_1 = center[1] + ((length / 2) * np.sin((cAngle)))
        endx_1 = center[0] + ((length / 2) * np.cos((cAngle)))
        endy_2 = center[1] - ((length / 2) * np.sin((cAngle)))
        endx_2 = center[0] - ((length / 2) * np.cos((cAngle)))
        start = np.asarray((endx_1, endy_1))
        end = np.asarray((endx_2, endy_2))
        max_direction, max_length = max(enumerate(abs(end - start)),
                                                  key=operator.itemgetter(1))
        dv = (end - start) / max_length
        line = [dda_round(start)]
        for step in range(int(max_length)):
            line.append(dda_round((step + 1) * dv + start))
        for coord in line:
            x_coord.append(coord[0])
            y_coord.append(coord[1])
            ts.append(i)
            pol.append(1)
    events = np.zeros((4, len(x_coord)))
    events[0, :] = np.asarray(x_coord)
    events[1, :] = np.asarray(y_coord)
    events[2, :] = np.asarray(ts)  
    events[3, :] = np.asarray(pol)

    ind = xy2ind(events[0, :], events[1, :], input_size_x)
    stimuli_generator = SpikeGeneratorGroup(input_size_x*input_size_y, indices=ind, 
                                    times=ts*ts_offset, name='rotating_bar')
    return stimuli_generator, events
