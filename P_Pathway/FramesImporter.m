
%% Importing Frames on Matlab
FramesFolder = 'C:\Users\RAX17\Desktop\Cartella condivisa Tesi\OUT spiking\Frames';
Dataset = '\Moving Bar';
xRes=240;
yRes=180;
FramesL={};
FramesR={};
xInputSize = 40;
yInputSize = 80;

count=0;
while(exist([FramesFolder Dataset '\R\' num2str(count) '.png'], 'file') == 2....
    && exist([FramesFolder Dataset '\L\' num2str(count) '.png'], 'file') == 2)
tmpFramesL=imread([FramesFolder Dataset '\L\' num2str(count) '.png']);
tmpFramesR=imread([FramesFolder Dataset '\R\' num2str(count) '.png']);
%FramesL{count+1}=tmpFramesL((yRes/2)-yInputSize/2:(yRes/2)+yInputSize/2,...
%    (xRes/2)-xInputSize/2:(xRes/2)+xInputSize/2);
%FramesR{count+1}=tmpFramesR((yRes/2)-yInputSize/2:(yRes/2)+yInputSize/2,...
%    (xRes/2)-xInputSize/2:(xRes/2)+xInputSize/2);
FramesL{count+1}=tmpFramesL;
FramesR{count+1}=tmpFramesR;
count = count + 1;
end


%% Saving Data

save(['DATA\' Dataset(2:end) 'FramesR.mat'],'FramesR');
save(['DATA\' Dataset(2:end) 'FramesL.mat'],'FramesL');
