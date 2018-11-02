clear
close all
clc
tic
addpath FUNCTIONS
addpath DATA
%%  POPULATION INIT
n_scales = 5;           % NUMBER OF SCALE - SET ON IMAGE SIZE AND DISPARITY CONTENT: max scale number n_scale, so that (minimum image size) / (2^n_scale) > (2*filter_size)
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



%% Load vL vR

Dataset = '/Moving Bar';

load(['DATA/' Dataset(2:end) 'FramesR.mat']);
load(['DATA/' Dataset(2:end) 'FramesL.mat']);

Out=cell(length(FramesL),1);

%% DISPARITY COMPUTATION
for ii=1:length(FramesL)
display(['Elaborating couple of Frames number '  num2str(ii)]);

tic;
    % CONVERSION TO DOUBLE
    II{ii}=cat(3,double(FramesL{ii}),double(FramesR{ii}));
    
    % POPULATION CODE
%     I=round(rand([250,400+15]));
%     II=I(:,1:400);
%     II(:,:,2)=I(:,16:end);
    Out{ii} = ctf_pop_disparity(II{ii},n_scales,n_filters,energy_th,ori_thr,ph_shift_file,filter_file);
    
   % figure,imagesc(D(:,:,1));
   % set(gca,'clim',[-30 30]);
   % title(energy_th);
    %clear L R II D
TimeperFrame=toc;
display(['Elaboration ended in ' num2str(TimeperFrame) ' seconds']);
display(['Extimated time left ' num2str(TimeperFrame*(length(FramesL)-ii)) ' seconds']);
display(['Couple of Frames left '  num2str(length(FramesL)-ii)]);

end


save Out %% � il vettore delle disparit� ha 2 mappe per i valori x e y 
% permette quindi di ottenere un vettore risultante delle disparit�
