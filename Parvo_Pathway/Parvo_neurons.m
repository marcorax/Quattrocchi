%% Parvo neurons, it refines the disparity computation from Magno neurons 
% while producing a dense estimantion.

clear all
clc;
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
% defined below (The code is meant ti work only with even values):
xInputSize = 40;
yInputSize = 80;

Dataset = '/Moving_Bar';
FramesFolder = '../Data/Extracted_Frames';
SpikesFolder = '../Data/Parvo_Population_Data' ;
load([FramesFolder Dataset '/FramesL.mat']); %Loading FrameL matrices
load([FramesFolder Dataset '/FramesR.mat']); %Loading FrameR matrices
load([FramesFolder Dataset '/Frametime.mat']); %Loading timeL and timeR arrays
load([SpikesFolder '/Input_On_' Dataset(2:end) '.mat']);
load([SpikesFolder '/Input_Off_' Dataset(2:end) '.mat']);

OnPop = {OnPop0,OnPop1,OnPop2,OnPop3,OnPop4,OnPop5,OnPop6,OnPop7};
OnPopt = {OnPop0t,OnPop1t,OnPop2t,OnPop3t,OnPop4t,OnPop5t,OnPop6t,OnPop7t};
OffPop = {OffPop0,OffPop1,OffPop2,OffPop3,OffPop4,OffPop5,OffPop6,OffPop7};
OffPopt = {OffPop0t,OffPop1t,OffPop2t,OffPop3t,OffPop4t,OffPop5t,OffPop6t,OffPop7t};
timeL = double(timeL)/1000000;
timeR = double(timeR)/1000000;

% I concatenate the Frames and resize them to focus the computation only 
% within the borders defined with xInputSize and yInputSize
shrinkdFramesL = cat(3,FramesL{:});
shrinkdFramesR = cat(3,FramesR{:});
shrinkdFramesL = shrinkdFramesL((yRes/2)-(yInputSize/2)+1:(yRes/2)+(yInputSize/2),...
                                (xRes/2)-(xInputSize/2)+1:(xRes/2)+(xInputSize/2),:);
shrinkdFramesR = shrinkdFramesR((yRes/2)-(yInputSize/2)+1:(yRes/2)+(yInputSize/2),...
                                (xRes/2)-(xInputSize/2)+1:(xRes/2)+(xInputSize/2),:);

%% P_Pathway segmentation
tic;
segmentate(shrinkdFramesL, shrinkdFramesR, OnPop, OnPopt, OffPop, OffPopt, timeL, timeOffs, NOrient)
toc;
%% M_Surface 
% Computing the disparity surfaces with information coming from the m
% pathway.
% This solution is obtained decoding disparity from the magno cells
% activity, a solution that ease the work done in this path but is not biological
% plausible. In future works, a bioinspired technique should be implemented. 

tic;
spkThr=3;
borders='nearest';
%borders='none';
[F_h, F_v, ScatterPerimeter] = m_surface(OnPop, OnPopt, OffPop, OffPopt,...
                                         timeL, timeOffs, NOrient, spkThr, ...
                                         borders);
toc

%% Picking up a frame + visualization 
iframe = 93;                                    
figure('Name','Surfaces by M Path','NumberTitle','off')
imgL = FramesL{iframe};
imgR = FramesR{iframe};
[xq,yq]=meshgrid(1:xInputSize,1:yInputSize);
%%%Left Image
colormap(gray(65536)); image(imgL); hold on;
pos = [(xRes/2-xInputSize/2),(yRes/2-yInputSize/2),xInputSize,yInputSize];
%adding firing cells
plot(ScatterPerimeter{iframe}(:,1)+(xRes-xInputSize)/2,...
ScatterPerimeter{iframe}(:,2)+(yRes-yInputSize)/2,'mo')
rectangle('Position',pos,'EdgeColor','r')
title('Left Frame')
hold off

%%%Right Image
figure('Name','Surfaces by M Path','NumberTitle','off')
colormap(gray(65536)); image(imgR); hold on;
pos = [(xRes/2-xInputSize/2),(yRes/2-yInputSize/2),xInputSize,yInputSize];
%adding estimations
xEstim = F_h{iframe}(ScatterPerimeter{iframe}(:,1),ScatterPerimeter{iframe}(:,2));
yEstim = F_v{iframe}(ScatterPerimeter{iframe}(:,1),ScatterPerimeter{iframe}(:,2));
xEstim = xEstim + ScatterPerimeter{iframe}(:,1)+(xRes-xInputSize)/2;
yEstim = yEstim + ScatterPerimeter{iframe}(:,2)+(yRes-yInputSize)/2;
plot(xEstim,yEstim,'bo')
rectangle('Position',pos,'EdgeColor','r')
title('Right Frame')
hold off;

%%%Horizontal Disparity surface
figure('Name','Surfaces by M Path','NumberTitle','off')
plot3(ScatterPerimeter{iframe}(:,1),ScatterPerimeter{iframe}(:,2)...
    ,ScatterPerimeter{iframe}(:,3),'mo')
colormap('default')
colorbar();
hold on
h_surface = F_h{iframe}(xq,yq);
%h_surface = imfilter(h_surface,h);
mesh(xq,yq,h_surface)
ax = gca;
ax.YDir = 'reverse';
axis([-20 60 0 80 -20 20])
title('M Surface Stream(Horizontal disparity)')
%legend('Sample Points','Interpolated Surface','Location','NorthWest')
hold off

%%%Vertical Disparity surface
figure('Name','Surfaces by M Path','NumberTitle','off')
plot3(ScatterPerimeter{iframe}(:,1),ScatterPerimeter{iframe}(:,2)...
    ,ScatterPerimeter{iframe}(:,4),'mo')
colormap('default')
colorbar;
hold on
v_surface = F_v{iframe}(xq,yq);
%v_surface = imfilter(v_surface,h);
mesh(xq,yq,F_v{iframe}(xq,yq))
ax = gca;
ax.YDir = 'reverse';
axis([-20 60 0 80 -20 20])
title('M Surface Stream(Vertical disparity)')
%legend('Sample Points','Interpolated Surface','Location','NorthWest')
hold off

x_disparity = F_h{iframe}(xq,yq);
y_disparity = F_v{iframe}(xq,yq);
%% Image morhping
x_disparity=medfilt2(x_disparity);
y_disparity=medfilt2(y_disparity);
 
x_disparity(isnan(x_disparity))=0;
y_disparity(isnan(y_disparity))=0;

% image filter
%  h = ones(5,5)/25;
%  x_disparity=imfilter(x_disparity,h);
%  y_disparity=imfilter(y_disparity,h);
 


imgL=imgL(((yRes/2)-yInputSize/2)+1:((yRes/2)+yInputSize/2),...
    ((xRes/2)-xInputSize/2)+1:((xRes/2)+xInputSize/2)); 
imgR=imgR(((yRes/2)-yInputSize/2)+1:((yRes/2)+yInputSize/2),...
    ((xRes/2)-xInputSize/2)+1:((xRes/2)+xInputSize/2));

OriginalImgL=imgL;
OriginalImgR=imgR;

mrfimgL=OriginalImgL;

for k = 1:length(xq(:))
    tmpxIndex = int8(max(min(xq(k)+x_disparity(yq(k),xq(k)),40),1));
    tmpyIndex = int8(max(min(yq(k)+y_disparity(yq(k),xq(k)),80),1));
    mrfimgL(yq(k),xq(k))=OriginalImgR(tmpyIndex,tmpxIndex);
end

for k = 1:length(xq(:))
    tmpxIndex = int8(max(min(xq(k)-x_disparity(yq(k),xq(k)),40),1));
    tmpyIndex = int8(max(min(yq(k)-y_disparity(yq(k),xq(k)),80),1));
    imgL(yq(k),xq(k))=mrfimgL(tmpyIndex,tmpxIndex);
end

% % image filter
% h = ones(5,5)/25;
% imgL=imfilter(imgL,h);
% imgR=imfilter(imgR,h);

CoarseMorphImgL=imgL;
CoarseMorphImgR=imgR;



%% DISPARITY COMPUTATION
%  POPULATION INIT
n_scales = 1;           % NUMBER OF SCALE - SET ON IMAGE SIZE AND DISPARITY CONTENT: max scale number n_scale, so that (minimum image size) / (2^n_scale) > (2*filter_size)
                        %                                                            max disparity decoded = +/- [2^n_scales /(4*0.063)]
energy_th = 2.2e-5;     % ENERGY THRESHOLD - SET ON DISPARITY RESULT
n_filters = 8;          % FIXED
ori_thr = 0;            % FIXED

% FILTER 11x11
ph_shift_file='FILTERS/ph_shift_PH7_F0.25.mat';
filter_file='FILTERS/Gt11B0.0833f0.25.mat';
% % FILTER 43x43
% filter='FILTERS/Gt43B0.0208f0.063.mat';
% ph_shift='FILTERS/ph_shift_PH7_F0.063.mat';


tic;
    % CONVERSION TO DOUBLE
    II=cat(3,double(imgL),double(imgR));
    
    % POPULATION CODE
%     I=round(rand([250,400+15]));
%     II=I(:,1:400);
%     II(:,:,2)=I(:,16:end);
    Out = ctf_pop_disparity(II,n_scales,n_filters,energy_th,ori_thr,ph_shift_file,filter_file);
    
%     figure,imagesc(Out(:,:,1)); hold on;
%     set(gca,'clim',[-20 20]);
%     title(energy_th);
%     colorbar;
%     %rectangle('Position',pos,'EdgeColor','r')
%     hold off;

Out((isnan(Out)))=0;

%save Out %% � il vettore delle disparit� ha 2 mappe per i valori x e y 
% permette quindi di ottenere un vettore risultante delle disparit�
toc
%% Result of the combination
h = ones(3,3)/9;
Out(:,:,1) = imfilter(Out(:,:,1),h);
Out(:,:,2) = imfilter(Out(:,:,2),h);

x_disparity = x_disparity-Out(:,:,1);
y_disparity = y_disparity-Out(:,:,2);


Refined_hdisparity=F_h{iframe}(xq,yq)-Out(:,:,1);
Refined_vdisparity=F_v{iframe}(xq,yq)-Out(:,:,2);
mrfimgL=OriginalImgL;

for k = 1:length(xq(:))
    tmpxIndex = int8(max(min(xq(k)+x_disparity(yq(k),xq(k)),40),1));
    tmpyIndex = int8(max(min(yq(k)+y_disparity(yq(k),xq(k)),80),1));
    mrfimgL(yq(k),xq(k))=OriginalImgR(tmpyIndex,tmpxIndex);
end

for k = 1:length(xq(:))
    tmpxIndex = int8(max(min(xq(k)-x_disparity(yq(k),xq(k)),40),1));
    tmpyIndex = int8(max(min(yq(k)-y_disparity(yq(k),xq(k)),80),1));
    imgL(yq(k),xq(k))=mrfimgL(tmpyIndex,tmpxIndex);
end

%filtering
% h = ones(5,5)/25;
% imgL=imfilter(imgL,h);
% imgR=imfilter(imgR,h);



FineMorphImgL=imgL;
FineMorphImgR=imgR;


%% Final results plots


%%%Left Image
figure('Name','Surfaces extraction and refining','NumberTitle','off')
colormap(gray(65536));
image(FramesL{iframe}); hold on;
pos = [(xRes/2-xInputSize/2),(yRes/2-yInputSize/2),xInputSize,yInputSize];
%adding firing cells
plot(ScatterPerimeter{iframe}(:,1)+(xRes-xInputSize)/2,...
ScatterPerimeter{iframe}(:,2)+(yRes-yInputSize)/2,'mo')
rectangle('Position',pos,'EdgeColor','r')
title('Left Frame')
hold off

%%%Right Image
figure('Name','Surfaces extraction and refining','NumberTitle','off')
colormap(gray(65536)); image(FramesR{iframe}); hold on;
pos = [(xRes/2-xInputSize/2),(yRes/2-yInputSize/2),xInputSize,yInputSize];
%adding estimations
xEstim = F_h{iframe}(ScatterPerimeter{iframe}(:,1),ScatterPerimeter{iframe}(:,2));
yEstim = F_v{iframe}(ScatterPerimeter{iframe}(:,1),ScatterPerimeter{iframe}(:,2));
xEstim = xEstim + ScatterPerimeter{iframe}(:,1)+(xRes-xInputSize)/2;
yEstim = yEstim + ScatterPerimeter{iframe}(:,2)+(yRes-yInputSize)/2;
plot(xEstim,yEstim,'bo')
rectangle('Position',pos,'EdgeColor','r')
title('Right Frame')
hold off;

%%%Horizontal Coarse Disparity surface
figure('Name','Surfaces extraction and refining','NumberTitle','off')
plot3(ScatterPerimeter{iframe}(:,1),ScatterPerimeter{iframe}(:,2)...
    ,ScatterPerimeter{iframe}(:,3),'mo')
colormap('default')
colorbar();
hold on
[xq,yq]=meshgrid(1:xInputSize,1:yInputSize);
h_surface = F_h{iframe}(xq,yq);
%h_surface = imfilter(h_surface,h);
mesh(xq,yq, h_surface)
ax3.YDir = 'reverse';
axis([-20 60 0 80 -20 20])
title('M Surface Stream(Horizontal disparity)')
%legend('Sample Points','Interpolated Surface','Location','NorthWest')
hold off

%%%Vertical Coarse Disparity surface
figure('Name','Surfaces extraction and refining','NumberTitle','off')
plot3(ScatterPerimeter{iframe}(:,1),ScatterPerimeter{iframe}(:,2)...
    ,ScatterPerimeter{iframe}(:,4),'mo')
colormap('default')
colorbar;
hold on
[xq,yq]=meshgrid(1:xInputSize,1:yInputSize);
v_surface = F_v{iframe}(xq,yq);
%v_surface = imfilter(v_surface,h);
mesh(xq,yq, v_surface)
ax4.YDir = 'reverse';
axis([-20 60 0 80 -20 20])
title('M Surface Stream(Vertical disparity)')
%legend('Sample Points','Interpolated Surface','Location','NorthWest')
hold off


%%%Horizontal Refined Disparity surface
figure('Name','Surfaces extraction and refining','NumberTitle','off')
plot3(ScatterPerimeter{iframe}(:,1),ScatterPerimeter{iframe}(:,2)...
    ,ScatterPerimeter{iframe}(:,3),'mo')
colormap('default')
colorbar();
hold on
[xq,yq]=meshgrid(1:xInputSize,1:yInputSize);
mesh(xq,yq,Refined_hdisparity)
ax5.YDir = 'reverse';
axis([-20 60 0 80 -20 20])
title('Refined (Horizontal disparity)')
%legend('Sample Points','Interpolated Surface','Location','NorthWest')
hold off

%%%Vertical Refined Disparity surface
figure('Name','Surfaces extraction and refining','NumberTitle','off')
plot3(ScatterPerimeter{iframe}(:,1),ScatterPerimeter{iframe}(:,2)...
    ,ScatterPerimeter{iframe}(:,4),'mo')
colormap('default')
colorbar;
hold on
[xq,yq]=meshgrid(1:xInputSize,1:yInputSize);
mesh(xq,yq,Refined_vdisparity)
ax6.YDir = 'reverse';
axis([-20 60 0 80 -20 20])
title('Refined Stream(Vertical disparity)')
%legend('Sample Points','Interpolated Surface','Location','NorthWest')
hold off




%%%Left Original Image
figure('Name','Morphing results','NumberTitle','off')
colormap(gray(65536)); image(OriginalImgL); hold on;
title('Left original image'); hold off;

%%%Right Original Image
figure('Name','Morphing results','NumberTitle','off')
colormap(gray(65536)); image(OriginalImgR); hold on;
title('Right original image'); hold off;

%%%Difference
figure('Name','Morphing results','NumberTitle','off')
colormap('default'); mesh(xq,yq,(double(OriginalImgR)-double(OriginalImgL))/2^16); hold on;
colorbar; caxis([-0.5,0.5]);
title('Difference between frames'); hold off;
corrMax=sum(sum(double(OriginalImgR).*double(OriginalImgR)));
OriginalCorr = sum(sum(double(OriginalImgL(20:60,12:30)).*double(OriginalImgR(20:60,12:30))))/corrMax
OriginalDiff = sum(sum(abs(double(OriginalImgR(20:60,12:30))-double(OriginalImgL(20:60,12:30)))))

%%%Left coarse-morphed Image
figure('Name','Morphing results','NumberTitle','off')
colormap(gray(65536)); image(CoarseMorphImgL); hold on;
title('Left coarse-morphed image'); hold off;

%%%Difference
figure('Name','Morphing results','NumberTitle','off')
colormap('default'); mesh(xq,yq,(double(CoarseMorphImgR)-double(CoarseMorphImgL))/2^16); hold on;
colorbar; caxis([-0.5,0.5]);
title('Difference between frames'); hold off;
CoarseCorr = sum(sum(double(CoarseMorphImgL(20:60,12:30)).*double(CoarseMorphImgR(20:60,12:30))))/corrMax
CoarseDiff = sum(sum(abs(double(CoarseMorphImgR(20:60,12:30))-double(CoarseMorphImgL(20:60,12:30)))))

%%%Left fine-morphed Image
figure('Name','Morphing results','NumberTitle','off')
colormap(gray(65536)); image(FineMorphImgL); hold on;
title('Left fine-morphed image'); hold off;

%%%Difference
figure('Name','Morphing results','NumberTitle','off')
colormap('default'); mesh(xq,yq,(double(FineMorphImgR)-double(FineMorphImgL))/2^16); hold on;
colorbar;caxis([-0.5,0.5]);
title('Difference between frame'); hold off;
FineCorr = sum(sum(double(FineMorphImgL(20:60,12:30)).*double(FineMorphImgR(20:60,12:30))))/corrMax
FineDiff = sum(sum(abs(double(FineMorphImgR(20:60,12:30))-double(FineMorphImgL(20:60,12:30)))))




