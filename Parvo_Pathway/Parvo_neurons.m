%% Parvo neurons, it refines the disparity computation from Magno neurons 
% while producing a dense estimantion.

clear all
clc
addpath Functions
addpath Filters

%Settings and Data

NOrient=8; % number of oriented cells, FIXED the use of pre computed linear
           % filters to compute the gabor parvo cells.
           % Note that this number CANNOT differ from the same paramater set 
           % for the Magno simple neurons, therefore before computing the M
           % pathway check if it is = 8, in future works this limitation might
           % be removed, thus the code in the magno path is written to accept
           % any number of orientations while being fixed on 8 
           
timeOffs=0.001; % The two pathways need to exchange information when they 
                % are temporally sinchronised, thus i compute disparity 
                % from the Magno_Pathway using events coming from a time
                % window long 2*timeOffs centered at each frame timestamps

% The resolution of the sensor
xRes = 240;
yRes = 180;

% The Magno pathway simulation on Brian2 requires a lot of time, therefore
% only a small centered portion of the image is computed, which borders are 
% defined below (The code is meant to work only with even values):
xInputSize = 40;
yInputSize = 80;

Dataset = '/Moving_Bar'; % Whenever you want to perform computation on a new 
                         % dataset, name it here and it should be fine 
                         
FramesFolder = '../Data/Extracted_Frames';
SpikesFolder = '../Data/Parvo_Population_Data' ;
ParvoResultsFolder = '../Data/Parvo_Population_Data';
% Loading Frame arrays and timestamps arrays
load([FramesFolder Dataset '/FramesL.mat'], 'FramesL'); %Loading FrameL matrices
load([FramesFolder Dataset '/FramesR.mat'], 'FramesR'); %Loading FrameR matrices
load([FramesFolder Dataset '/Frametime.mat'], 'timeL', 'timeR'); %Loading timeL and timeR arrays

% Loading Magno disparity neuron activity, computed previously
load([SpikesFolder '/Input_On_' Dataset(2:end) '.mat'],'OnPop0', 'OnPop1',...
                'OnPop2', 'OnPop3', 'OnPop4', 'OnPop5', 'OnPop6', 'OnPop7',...
                'OnPop0t', 'OnPop1t', 'OnPop2t', 'OnPop3t', 'OnPop4t', 'OnPop5t',...
                'OnPop6t', 'OnPop7t');
load([SpikesFolder '/Input_Off_' Dataset(2:end) '.mat'],'OffPop0', 'OffPop1',...
                'OffPop2', 'OffPop3', 'OffPop4', 'OffPop5', 'OffPop6', 'OffPop7',...
                'OffPop0t', 'OffPop1t', 'OffPop2t', 'OffPop3t', 'OffPop4t', 'OffPop5t',...
                'OffPop6t', 'OffPop7t');

% On and off populations rearranging
OnPop = {OnPop0,OnPop1,OnPop2,OnPop3,OnPop4,OnPop5,OnPop6,OnPop7};
clear OnPop0 OnPop1 OnPop2 OnPop3 OnPop4 OnPop5 OnPop6 OnPop7
OnPopt = {OnPop0t,OnPop1t,OnPop2t,OnPop3t,OnPop4t,OnPop5t,OnPop6t,OnPop7t};
clear OnPop0t OnPop1t OnPop2t OnPop3t OnPop4t OnPop5t OnPop6t OnPop7t
OffPop = {OffPop0,OffPop1,OffPop2,OffPop3,OffPop4,OffPop5,OffPop6,OffPop7};
clear OffPop0 OffPop1 OffPop2 OffPop3 OffPop4 OffPop5 OffPop6 OffPop7
OffPopt = {OffPop0t,OffPop1t,OffPop2t,OffPop3t,OffPop4t,OffPop5t,OffPop6t,OffPop7t};
clear OffPop0t OffPop1t OffPop2t OffPop3t OffPop4t OffPop5t OffPop6t OffPop7t
timeL = double(timeL)/1000000;
timeR = double(timeR)/1000000;

% I concatenate the Frames and resize them to focus the computation only 
% within the borders defined with xInputSize and yInputSize
cutFramesL = cat(3,FramesL{:});
cutFramesR = cat(3,FramesR{:});
cutFramesL = cutFramesL((yRes/2)-(yInputSize/2)+1:(yRes/2)+(yInputSize/2),...
                                (xRes/2)-(xInputSize/2)+1:(xRes/2)+(xInputSize/2),:);
cutFramesR = cutFramesR((yRes/2)-(yInputSize/2)+1:(yRes/2)+(yInputSize/2),...
                                (xRes/2)-(xInputSize/2)+1:(xRes/2)+(xInputSize/2),:);

%% P_Pathway segmentation
% To avoid errors of interpolation performed by the segment below, the
% surfaces are segmentred with the informations extracted from the frames

% Number of iteration of the filling algorithm.
nIterations = 50;

% The border thikness (the size of each point generated by the active neurons
% in the mask used by activecontour)
bord_thick = 2;

% Note that the function will signal with a warning whenever the segmentation
% will not be possible (disabling warnings won't affect it beacuse we all know
% that MATLAB is an high level piece of shit ¯\_(ツ)_/¯)

tic;
[seg_frames_l,seg_frames_r] = segmentate(cutFramesL, cutFramesR, OnPop, OnPopt,...
                            OffPop, OffPopt, timeL, timeOffs, NOrient, ...
                            nIterations, bord_thick);
toc
%% M_Surface 
% Computing the disparity surfaces with information coming from the m_pathway.
% This solution is obtained decoding disparity from the magno cells
% activity, a solution that ease the work done in this path but is not biological
% plausible. In future works, a bioinspired technique should be implemented. 
% The surfaces are then segmented using the results from the P_Pathawy
% segmentation, only in respect to the left frame, if anything close to
% this computation is happening in the brain I find difficult to believe, 
% that the computation might happen unilateral as here. But probably
% informations are merged and reinforced building disparity surfaces
% respective for both eyes.

tic;
spkThr=3;
borders='linear';
[Coarse_h_disparity, Coarse_v_disparity, ScatterPerimeter] = m_surface(OnPop,...
                                        OnPopt, OffPop, OffPopt, timeL,...
                                        timeOffs, NOrient, spkThr, ...
                                        borders, seg_frames_l);
toc
%% Parvo population disparity computation and refining
tic;
% Polpulation settings
energy_th = 2.2e-5;     % ENERGY THRESHOLD - SET ON DISPARITY RESULT
n_filters = 8;          % FIXED (To remember that this number cannot be set)
ori_thr = 0;            % FIXED

% Receptive field sizes (Select by commenting the unwanted one)
% FILTER 11x11 
ph_shift_file='FILTERS/ph_shift_PH7_F0.25.mat';
filter_file='FILTERS/Gt11B0.0833f0.25.mat';
% % FILTER 43x43
% filter='FILTERS/Gt43B0.0208f0.063.mat';
% ph_shift='FILTERS/ph_shift_PH7_F0.063.mat';

% Preparing data for P_Disparity function
II=cat(4,double(cutFramesL),double(cutFramesR));
mask=cat(4,seg_frames_l,seg_frames_r);
coarse_disp=cat(4,Coarse_h_disparity,Coarse_v_disparity);

% Computing refined disparity
fine_disp = P_Disparity(II,coarse_disp,mask,energy_th,ori_thr,ph_shift_file,filter_file);    
toc
%% Save Results
save([ParvoResultsFolder '/Disparity_results_' Dataset(2:end) '.mat'], 'mask', 'coarse_disp', 'fine_disp', 'ScatterPerimeter')
%% Load Results 
load([ParvoResultsFolder '/Disparity_results_' Dataset(2:end) '.mat'], 'mask', 'coarse_disp', 'fine_disp', 'ScatterPerimeter')
%% Frame position for plotting
% Set the reference frame for the plots below
% Remember to load also the frames and frametime in the first subsection 
% of the script if you want to plot everything
frame = 93;

%% Reference frames computing window plot 
% The window where the computation has been performed (has set in inputSize values)

imgL = FramesL{frame};
imgR = FramesR{frame};
% Computing window
pos = [(xRes/2-xInputSize/2),(yRes/2-yInputSize/2),xInputSize,yInputSize];
%%% Left Image
figure('Name','Input frames','NumberTitle','off')
colormap(gray(65536)); image(imgL); hold on;
rectangle('Position',pos,'EdgeColor','r')
title('Left Frame')
hold off
%%% Left Image
figure('Name','Input frames','NumberTitle','off')
colormap(gray(65536)); image(imgR); hold on;
rectangle('Position',pos,'EdgeColor','r')
title('Right Frame')
hold off


%% Disparity plots

%%% Horizontal Coarse Disparity surface
figure('Name','Surfaces extraction and refining','NumberTitle','off')
% Left Magno Receptive Fields Positions (True values, not interpolated)
plot3(ScatterPerimeter{frame}(:,1),ScatterPerimeter{frame}(:,2)...
    ,ScatterPerimeter{frame}(:,3),'mo')
colormap('default')
colorbar();
hold on
[xq,yq]=meshgrid(1:xInputSize,1:yInputSize);
mesh(xq, yq, coarse_disp(:,:,frame,1))
ax5.YDir = 'reverse';
axis([-20 60 0 80 -20 20])
title('Coarse (Horizontal disparity)')
legend('Sample Points','Interpolated Surface','Location','NorthWest')
hold off

%%% Vertical Coarse Disparity surface
figure('Name','Surfaces extraction and refining','NumberTitle','off')
% Left Magno Receptive Fields Positions (True values, not interpolated)
plot3(ScatterPerimeter{frame}(:,1),ScatterPerimeter{frame}(:,2)...
    ,ScatterPerimeter{frame}(:,4),'mo')
colormap('default')
colorbar;
hold on
[xq,yq]=meshgrid(1:xInputSize,1:yInputSize);
mesh(xq, yq, coarse_disp(:,:,frame,2))
ax6.YDir = 'reverse';
axis([-20 60 0 80 -20 20])
title('Coarse (Vertical disparity)')
legend('Sample Points','Interpolated Surface','Location','NorthWest')
hold off

%%% Horizontal Coarse Disparity surface
figure('Name','Surfaces extraction and refining','NumberTitle','off')
% Left Magno Receptive Fields Positions (True values, not interpolated)
plot3(ScatterPerimeter{frame}(:,1),ScatterPerimeter{frame}(:,2)...
    ,ScatterPerimeter{frame}(:,3),'mo')
colormap('default')
colorbar();
hold on
[xq,yq]=meshgrid(1:xInputSize,1:yInputSize);
mesh(xq, yq, coarse_disp(:,:,frame,1)+fine_disp(:,:,frame,1))
ax5.YDir = 'reverse';
axis([-20 60 0 80 -20 20])
title('Fine (Horizontal disparity)')
legend('Sample Points','Interpolated Surface','Location','NorthWest')
hold off

%%% Vertical Coarse Disparity surface
figure('Name','Surfaces extraction and refining','NumberTitle','off')
% Left Magno Receptive Fields Positions (True values, not interpolated)
plot3(ScatterPerimeter{frame}(:,1),ScatterPerimeter{frame}(:,2)...
    ,ScatterPerimeter{frame}(:,4),'mo')
colormap('default')
colorbar;
hold on
[xq,yq]=meshgrid(1:xInputSize,1:yInputSize);
mesh(xq, yq, coarse_disp(:,:,frame,2)+fine_disp(:,:,frame,2))
ax6.YDir = 'reverse';
axis([-20 60 0 80 -20 20])
title('Fine (Vertical disparity)')
legend('Sample Points','Interpolated Surface','Location','NorthWest')
hold off

%% Morph plots (A qualitative test)
% These are the results of the morphing test,
% where synthetic left frames are built with original right frames pixel
% whose positions are computed using the disparity map previosly computed

% Original Right Frame
original_left_frame = cutFramesL(:,:,frame);
right_frame = cutFramesR(:,:,frame);
h_coarse_disp = squeeze(coarse_disp(:,:,frame,1));
v_coarse_disp = squeeze(coarse_disp(:,:,frame,2));
figure('Name','Morphing Test','NumberTitle','off')
colormap(gray(65536)); image(original_left_frame); hold on;
title('Left Frame (Original)')
hold off;

% Computing the pixel positions (coarse)
coarse_morphed_frame = original_left_frame;
[ly,lx]=size(original_left_frame);
[indy,indx] = find(abs(h_coarse_disp)+abs(v_coarse_disp)~=0);
indy_shift = int32(indy+round(v_coarse_disp(sub2ind([ly,lx],indy,indx))));
indx_shift = int32(indx+round(h_coarse_disp(sub2ind([ly,lx],indy,indx))));
% remove all shifts outside boundaries
valid_ind = (0<indy_shift)&(indy_shift<=ly)&(0<indx_shift)&(indx_shift<=lx);
indy_shift = indy_shift(valid_ind);
indx_shift = indx_shift(valid_ind);
indy = indy(valid_ind);
indx = indx(valid_ind);
coarse_morphed_frame(sub2ind([ly,lx],indy,indx)) = ...
    right_frame(sub2ind([ly,lx],indy_shift,indx_shift));

% Morphed Right Frame (coarse)
figure('Name','Morphing Test','NumberTitle','off')
colormap(gray(65536)); image(coarse_morphed_frame); hold on;
title('Morphed Left Frame (coarse)')
hold off;

% Computing the pixel positions (fine)
fine_morphed_frame = original_left_frame;
refined_disp = coarse_disp+fine_disp;
h_fine_disp = squeeze(refined_disp(:,:,frame,1));
v_fine_disp = squeeze(refined_disp(:,:,frame,2));
[indy,indx] = find(abs(h_fine_disp)+abs(v_fine_disp)~=0);
indy_shift = int32(indy+round(v_fine_disp(sub2ind([ly,lx],indy,indx))));
indx_shift = int32(indx+round(h_fine_disp(sub2ind([ly,lx],indy,indx))));
% remove all shifts outside boundaries
valid_ind = (0<indy_shift)&(indy_shift<=ly)&(0<indx_shift)&(indx_shift<=lx);
indy_shift = indy_shift(valid_ind);
indx_shift = indx_shift(valid_ind);
indy = indy(valid_ind);
indx = indx(valid_ind);
fine_morphed_frame(sub2ind([ly,lx],indy,indx)) = ...
    right_frame(sub2ind([ly,lx],indy_shift,indx_shift));

% Morphed Right Frame (fine)
figure('Name','Morphing Test','NumberTitle','off')
colormap(gray(65536)); image(fine_morphed_frame); hold on;
title('Morphed Left Frame (fine)')
hold off;

% Difference between coarse and original left frame
figure('Name','Morphing Test','NumberTitle','off')
colormap('default'); mesh(xq,yq,(double(original_left_frame)-double(coarse_morphed_frame))/2^16); hold on;
colorbar; caxis([-0.5,0.5]);
title('Difference between coarse morphed and original image'); hold off;

corrMax=sum(sum(double(coarse_morphed_frame).*double(coarse_morphed_frame)));
'Correlation between coarse morphed and original image'
coarseCorr = sum(sum(double(original_left_frame).*double(coarse_morphed_frame)))/corrMax

% Difference between coarse and original left frame
figure('Name','Morphing Test','NumberTitle','off')
colormap('default'); mesh(xq,yq,(double(original_left_frame)-double(fine_morphed_frame))/2^16); hold on;
colorbar; caxis([-0.5,0.5]);
title('Difference between fine morphed and original image'); hold off;

corrMax=sum(sum(double(fine_morphed_frame).*double(fine_morphed_frame)));
'Correlation between fine morphed and original image'
fineCorr = sum(sum(double(original_left_frame).*double(fine_morphed_frame)))/corrMax
