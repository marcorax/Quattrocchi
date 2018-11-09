# -*- coding: utf-8 -*-
# @Author: marcorax93

''' 
This file contains functions used to connect and set weights in Brian2, and tili
libraries (to know more, please take a look to the functions description in the 
Brian2 guide https://brian2.readthedocs.io/en/2.0rc/user/functions.html)

Note that kernels here are meant to be connected a spike generator, this raise 
the issue that instead of input neurons with x and y coordinates we need to account 
for linear i coordinates.
I hope that this issue will be soon solved, and the spike generator will allow
to address the spikes with x and y like a normal neuron group.
(in that case the functions won't need the i index anympore as input but they 
will use directly x_pre and y_pre as the post synaptic counterparts)

The positional values (x_post, y_post) are all expected to be 
indeces of a 2D network (positive integer)
These function are not expected to work with complex and neuro-topological inspired 
models

'''
import numpy as np
from brian2 import implementation, check_units, ms, exp, mean, diff, declare_types


@implementation('cpp', '''
int connect_spike_gen(int i, int x_post, int y_post, int input_size_x, int input_size_y, int kernel_size){
    int ix = i % input_size_x;
    int iy = i / input_size_x;
    int result = (abs(ix - x_post)<kernel_size/2)*(abs(iy - y_post)<kernel_size/2);
    return result;
    }
     ''')
@declare_types(i='integer', x_post='integer', y_post='integer', input_size_x='integer', input_size_y='integer', kernel_size='integer', result='integer')
@check_units(i=1, x_post=1, y_post=1, input_size_x=1, input_size_y=1, kernel_size=1, result=1)
def connect_spike_gen(i, x_post, y_post, input_size_x, input_size_y, kernel_size):
    """Summary: function used to connect a cortical population with a spatial
    rapresentation (x and y indices) and a (spike generator)
    with a squared kernel of connections of size kernel_size**2
    

    Args:
        i (int): Spike generator index
        x_post (int): x index of the cortical cell
        y_post (int): y index of the cortical cell
        input_size_x: horizontal dimension of the retinal plane
        input_size_y: vertical dimension of the retinal plane
        kernel_size (int): linear dimension of the squared connectivity kernel,
                            it must be an odd number

    Returns:
        int: 0 => the selected neurons are not connected, 
             1 => the selected neurons are connecteds
    """

    (iy, ix) = np.unravel_index(i, (input_size_y, input_size_x))
    result = (abs(ix - x_post)<kernel_size/2) *(abs(iy - y_post)<kernel_size/2)
    return result




@implementation('cpp', '''
float gabor_spike_gen(int i, int x_post, int y_post, int input_size_x, int input_size_y, float theta, float sigma_x, float sigma_y, float freq){
    int ix = i % input_size_x;
    int iy = i / input_size_x;
    float x =  (ix - x_post)*cos(theta-M_PI/2) + (iy - y_post)*sin(theta-M_PI/2);
    float y = -(ix - x_post)*sin(theta-M_PI/2) + (iy - y_post)*cos(theta-M_PI/2);
    float exponent = -((pow(x,2)/(2*pow(sigma_x,2))) + (pow(y,1)/(2*pow(sigma_y,2))));
    float result = exp(exponent)*cos(M_PI*x/freq);
    return result;
    }
    
     ''')
@declare_types(i='integer', x_post='integer', y_post='integer', input_size_x='integer', input_size_y='integer', theta='float', sigma_x='float', sigma_y='float', freq='float', result='float')
@check_units(i=1, x_post=1, y_post=1, input_size_x=1, input_size_y=1, theta=1, sigma_x=1, sigma_y=1, freq=1, result=1)
def gabor_spike_gen(i, x_post, y_post, input_size_x, input_size_y,  theta, sigma_x, sigma_y, freq):
    """Summary: function used to compute the weights of simple cells connected 
    to an artificial  (the spike generator), please note that this gabor 
    rapresentation is not indicative of any actual connectivity pattern in the
    brain because it doesn't account for multiple layers, retinocortical 
    trasformations, or exclusivity costraints on excitations or hinibitions
    (excitatory-hinibitory only cells )
    

    Args:
        i (int): Spike generator index
        x_post (int): x index of the cortical cell
        y_post (int): y index of the cortical cell
        input_size_x: horizontal dimension of the retinal plane
        input_size_y: vertical dimension of the retinal plane
        theta (float): orientation of the gabor connectivity kernel
        sigma_x (float): variance along x 
        sigma:y (float): variance along y
        freq (float): frequency of the filter

    Returns:
        float: The weight between the i and j neuron

    """

    (iy, ix) = np.unravel_index(i, (input_size_y, input_size_x))
    x =  (ix - x_post)*np.cos(theta-np.pi/2) + (iy - y_post)*np.sin(theta-np.pi/2)
    y = -(ix - x_post)*np.sin(theta-np.pi/2) + (iy - y_post)*np.cos(theta-np.pi/2)
    exponent = -(((x**2)/(2*sigma_x**2)) + ((y**2)/(2*sigma_y**2)))
    result = exp(exponent)*np.cos(np.pi*x/freq)
    return result


@implementation('cpp', '''
float gaussian_spike_gen(int i, int x_post, int y_post, int input_size_x, int input_size_y, float sigma){
    int ix = i % input_size_x;
    int iy = i / input_size_x;    
    float x = ix - x_post;
    float y = iy - y_post;
    float exponent = -((pow(x,2)+pow(y,2))/(2*pow(sigma,2)));
    float result = exp(exponent);
    return result;
    }
    
     ''')
@declare_types(i='integer', x_post='integer', y_post='integer', input_size_x='integer', input_size_y='integer', sigma='float', result='float')
@check_units(i=1, x_post=1, y_post=1, input_size_x=1, input_size_y=1, sigma=1, result=1)
def gaussian_spike_gen(i, x_post, y_post, input_size_x, input_size_y, sigma):
    """Summary: function used to compute the weights of synapses of neurons 
        connected to a linear spike generator group in a 2D gaussian configuration 

    Args:
        i (int): Spike generator index
        x_post (int): x index of the post_synaptic cell
        y_post (int): y index of the post_synaptic cell
        input_size_x: horizontal dimension of the retinal plane
        input_size_y: vertical dimension of the retinal plane
        sigma (float): variance along the radius


    Returns:
        float: The weight between the i and j neuron

    """

    (iy, ix) = np.unravel_index(i, (input_size_y, input_size_x))
    x = ix - x_post
    y = iy - y_post
    exponent = -(((x**2)+(y**2))/(2*sigma**2))
    result = exp(exponent)
    return result


def stereo_centers(input_size_x, input_size_y, spacing, kernel_size, depth_planes):
    """Summary: function used to compute the connectivity kernel centers for 
    stereo cells, used for setting xl, xr, yl, yr of stereo cells in the net

    Args:
        
        kernel_size (int) : The linear size of the squared connectivity kernels
        spacing (int) : The number of pixels between two consequent kernel centers,
                        if 0, the first RF will be in (a,a) and the second will
                        be positioned in (a+1,a)
        depth_planes (int) : The number of disparity levels which the stereo cells
                            encode, must be an odd number (the zero disparity 
                            level must to be taken in account)
    
    Returns:
        neuron_centers [[xl],[yl],[xr],[yr]] (int) : list including all the 
                                                    connectivity cell centers
        num_stereo_cells (int) : number of the total amount of stereo cells recquired
                                 to build the network
        """
    num_ker_x = int(np.ceil((input_size_x-(kernel_size-1))/(1+spacing)))
    num_ker_y = int(np.ceil((input_size_y-(kernel_size-1))/(1+spacing)))
    neuron_centers = np.zeros((num_ker_x*num_ker_y*depth_planes*depth_planes,4),int)
    tmpY,tmpX = np.dot(np.unravel_index(range(depth_planes**2),(depth_planes,depth_planes)),(spacing+1))
    tmpX=tmpX-int(depth_planes/2)*(spacing+1)
    tmpY=tmpY-int(depth_planes/2)*(spacing+1)

    for i in range(num_ker_y):
        for j in range(num_ker_x):
            lindex=(j+i*num_ker_x)*depth_planes**2
            rindex=lindex+depth_planes**2
            centerX=int((j*(spacing+1))+int(kernel_size-1)/2)
            centerY=int((i*(spacing+1))+int(kernel_size-1)/2)
            neuron_centers[lindex:rindex,0:2]=[centerX,centerY]
            neuron_centers[lindex:rindex,2:4] = np.transpose([ tmpX + centerX, tmpY + centerY])

    #I want to remove every stereocell with RF Outside the Image Space
    
    #boolean array where i identify x values included in the InputSizeX (considering the borders which cannot 
    #contain another RF)
    
    full_kernels_x=(((kernel_size-1)/2)<=neuron_centers[:,[0,2]])*(neuron_centers[:,[0,2]]<(input_size_x-((kernel_size-1)/2)))
    full_kernels_y=(((kernel_size-1)/2)<=neuron_centers[:,[1,3]])*(neuron_centers[:,[1,3]]<(input_size_y-((kernel_size-1)/2)))
    
    neuron_centers=neuron_centers[np.all(full_kernels_x*full_kernels_y,axis=1),:]
    num_stereo_cells=np.shape(neuron_centers)[0]
    return(neuron_centers, num_stereo_cells)