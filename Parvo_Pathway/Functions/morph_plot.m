%% Picking up a frame + visualization 
iframe = 93;                                    
figure('Name','Surfaces by M Path','NumberTitle','off')
imgL = FramesL{iframe};
imgR = FramesR{iframe};
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
% Estimations, left disparity neurons receptive fields centers + computed 
% disparity
Index=sub2ind(size(Coarse_h_disparity{iframe}),ScatterPerimeter{iframe}(:,2),ScatterPerimeter{iframe}(:,1));
xEstim = Coarse_h_disparity{iframe}(Index);
yEstim = Coarse_v_disparity{iframe}(Index);
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
mesh(Coarse_h_disparity{iframe})
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
mesh(Coarse_v_disparity{iframe})
ax = gca;
ax.YDir = 'reverse';
axis([-20 60 0 80 -20 20])
title('M Surface Stream(Vertical disparity)')
%legend('Sample Points','Interpolated Surface','Location','NorthWest')
hold off

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

