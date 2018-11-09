# -*- coding: utf-8 -*-
# @Author: marcorax93

''' 
This script decodes aedat file and extracts frames,
that saves in png and .mat files

'''
import os
import struct
import time
import numpy as np
import png
from Davisloading import framesToPy
import scipy.io as sio


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return (idx)


aedatRecordR = "/home/marcorax93/Documenti/DAVIS Recordings/Frames_R_Moving_Bar-2018_03_06_17_04_05.aedat" 
outfolderR = "/home/marcorax93/Documenti/Disparitydetector Brian2/Frames/Moving Bar/R/"

aedatRecordL = "/home/marcorax93/Documenti/DAVIS Recordings/Frames_L_Moving_Bar-2018_03_06_17_04_05.aedat" 
outfolderL = "/home/marcorax93/Documenti/Disparitydetector Brian2/Frames/Moving Bar/L/"

xCameraRes=240
yCameraRes=180

aedatvidL=framesToPy(filename=aedatRecordL, xdim=xCameraRes, ydim=yCameraRes)
aedatvidR=framesToPy(filename=aedatRecordR, xdim=xCameraRes, ydim=yCameraRes)

framesR=[]
framesL=[]

timeR=[]
timeL=[]

if (len(aedatvidR[0])>len(aedatvidL[0])):
    for k in aedatvidL[2]:
        index=find_nearest(aedatvidR[2],k)
        framesR.append(aedatvidR[0][index])
        timeR.append(aedatvidR[2][index])
    framesL = aedatvidL[0]
    timeL = aedatvidL[2]
elif (len(aedatvidL[0])>len(aedatvidR[0])):
    for k in aedatvidR[2]:
        index=find_nearest(aedatvidL[2],k)
        framesL.append(aedatvidL[0][index])
        timeL.append(aedatvidL[2][index])
    framesR = aedatvidR[0]
    timeR = aedatvidR[2]
else:
    framesL = aedatvidL[0]
    framesR = aedatvidR[0]
    timeL = aedatvidL[2]
    timeR = aedatvidR[2]


#%% Saving 


framecounter=0
for index in np.argsort(timeL):
    with open(outfolderL+str(framecounter)+'.png', 'wb') as f:
        writer = png.Writer(width=xCameraRes, height=yCameraRes, bitdepth=16, greyscale=True)
        zgray2list = framesL[index]
        writer.write(f, zgray2list)
    framecounter += 1


framecounter=0
for index in np.argsort(timeR):
    with open(outfolderR+str(framecounter)+'.png', 'wb') as f:
        writer = png.Writer(width=xCameraRes, height=yCameraRes, bitdepth=16, greyscale=True)
        zgray2list = framesR[index]
        writer.write(f, zgray2list)
    framecounter += 1
    
        
sio.savemat('Frametime.mat', {'timeL':timeL,'timeR':timeR})   
