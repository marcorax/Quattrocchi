%% Energy model
% Originary energy model working only on frames, it is here implemented
% to perform testing with the Quattrocchi model.

clear all
clc
addpath Functions
addpath Filters



% The resolution of the sensor
xRes = 240;
yRes = 180;

% To compare apples to apples is advise to use the same computation window
% as used in Parvo Neurons
xInputSize = 40;
yInputSize = 80;

%Settings and Data

n_scales = 2;           % NUMBER OF SCALE - SET ON IMAGE SIZE AND DISPARITY CONTENT: max scale number n_scale, so that (minimum image size) / (2^n_scale) > (2*filter_size)
                        %                                                            max disparity decoded = +/- [2^n_scales /(4*0.063)]
energy_th = 2.2e-5;     % ENERGY THRESHOLD - SET ON DISPARITY RESULT
n_filters = 8;          % FIXED
ori_thr = 0;            % FIXED

% Receptive field sizes (Select by commenting the unwanted one)
% FILTER 11x11
ph_shift_file='Filters/ph_shift_PH7_F0.25.mat';
filter_file='Filters/Gt11B0.0833f0.25.mat';
% % FILTER 43x43
% filter='Filters/Gt43B0.0208f0.063.mat';
% ph_shift='Filters/ph_shift_PH7_F0.063.mat';

Dataset = '/Moving_Bar'; % Whenever you want to perform computation on a new 
                         % dataset, name it here and it should be fine 
                         
FramesFolder = '../Data/Extracted_Frames';
ParvoResultsFolder = '../Data/Parvo_Population_Data';
% Loading Frame arrays and timestamps arrays
load([FramesFolder Dataset '/FramesL.mat'], 'FramesL'); %Loading FrameL matrices
load([FramesFolder Dataset '/FramesR.mat'], 'FramesR'); %Loading FrameR matrices
load([FramesFolder Dataset '/Frametime.mat'], 'timeL', 'timeR'); %Loading timeL and timeR arrays

% I concatenate the Frames and resize them to focus the computation only 
% within the borders defined with xInputSize and yInputSize
cutFramesL = cat(3,FramesL{:});
cutFramesR = cat(3,FramesR{:});
cutFramesL = cutFramesL((yRes/2)-(yInputSize/2)+1:(yRes/2)+(yInputSize/2),...
                                (xRes/2)-(xInputSize/2)+1:(xRes/2)+(xInputSize/2),:);
cutFramesR = cutFramesR((yRes/2)-(yInputSize/2)+1:(yRes/2)+(yInputSize/2),...
                                (xRes/2)-(xInputSize/2)+1:(xRes/2)+(xInputSize/2),:);
%% Disparity computation
tic;
% Preparing data to be computed
II=cat(4,double(cutFramesL),double(cutFramesR));

% Computing disparity
disparity = ctf_pop_disparity(II,n_scales,energy_th,ori_thr,ph_shift_file,filter_file);    
toc

%% Save Results
save([ParvoResultsFolder '/Energy_model_Disparity_results_' Dataset(2:end) '.mat'], 'disparity')
%% Load Results 
load([ParvoResultsFolder '/Energy_model_Disparity_results_' Dataset(2:end) '.mat'], 'disparity')

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

%%% Horizontal Disparity 
figure('Name','Disparity results','NumberTitle','off')
colormap('default')
colorbar();
hold on
[xq,yq]=meshgrid(1:xInputSize,1:yInputSize);
mesh(xq, yq, disparity(:,:,frame,1))
ax5.YDir = 'reverse';
axis([-20 60 0 80 -20 20])
title('Horizontal disparity')
hold off

%%% Vertical Disparity
figure('Name','Disparity results','NumberTitle','off')
colormap('default')
colorbar;
hold on
[xq,yq]=meshgrid(1:xInputSize,1:yInputSize);
mesh(xq, yq, disparity(:,:,frame,2))
ax6.YDir = 'reverse';
axis([-20 60 0 80 -20 20])
title('Vertical disparity')
hold off


