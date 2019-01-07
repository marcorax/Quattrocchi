function Disparity = ctf_pop_disparity(frames,n_scales,energy_thr,ori_thr,ph_shift_file,filter_file)

[sy1, sx1, n_frames, st] = size(frames);
D = cell(length(n_frames),1);
parfor frame = 1:n_frames
    %%%%%%%%%%%%%%%%%%%%%%%
    % PADDING FOR PYRAMID %
    %%%%%%%%%%%%%%%%%%%%%%%
    G = cell(2);
    II = squeeze(frames(:,:,frame,:));
    
    if (st~=2)
        error('size(II,3) should equal 2!\n');
    end
    fac = 2^(n_scales-1);

    sy2 = ceil(sy1 ./ fac) .* fac; % target resolution
    sx2 = ceil(sx1 ./ fac) .* fac; % target resolution

    II = [ II ; repmat(II(end,:,:),[sy2-sy1 1 1]) ]; % replicate border row
    II = [ II repmat(II(:,end,:),[1 sx2-sx1 1]) ];   % replicate border column

    II = image_pyramid(II,n_scales);


    %%%%%%%%%%%
    % Level 1 %
    %%%%%%%%%%%
    % display('Scale #: 1');

    % COMPUTE GABOR FILTERING ON COARSER SCALE (REAL AND IMAG PART)
    [F] = filt_gabor(II{1},filter_file);

    % COMPUTE PHASE SHIFT ON COARSER SCALE
    G{1,1}=F{1}(:,:,:,1); G{1,2}=F{2}(:,:,:,1);
    G(2,:)=shift_in_phase(F{1}(:,:,:,2),F{2}(:,:,:,2),ph_shift_file);
    

    % PLOT FILTER SHAPE (uncomment to use it)
    % plot_filter(G)

    % POPULATION DECODING
    tmp_disp = population(G,energy_thr,ori_thr,ph_shift_file);

    % MASK UNRELIABLE VALUE
    invalid = isnan(tmp_disp);

    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Coarse-to-fine Estimation and Merging %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    for scale = 2:n_scales

        % display(['Scale #: ' num2str(scale)]);

        % EXPAND DISPARITY OF PREVIOUS SCALE
        tmp_disp = expand(tmp_disp.*2);    

        % COMPUTE GABOR FILTERING ON CURRENT SCALE (REAL AND IMAG PART)
        [F] = filt_gabor(II{scale},filter_file);

        % COMPUTE PHASE SHIFT ON CURRENT SCALE
        G{1,1}=F{1}(:,:,:,1); G{1,2}=F{2}(:,:,:,1);
        G(2,:)=shift_in_phase(F{1}(:,:,:,2),F{2}(:,:,:,2),ph_shift_file);
        

        % WARP POPULATION RESPONSE ACCORDING TO DISPARITY OF PREVIOUS SCALE
        G = warp_sequence(G,tmp_disp);

        % POPULATION DECODING
        Ds = population(G,energy_thr,ori_thr,ph_shift_file);

        % MERGE DISPARITY OF PREVIOUS AND CURRENT SCALE
        tmp_disp = merge_disparity(tmp_disp,Ds);

        % MASK UNRELIABLE VALUES
        invalid = isnan(Ds);

        
    end

    % MASK UNRELIABLE VALUE
    tmp_disp(invalid) = NaN;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REMOVE PADDING FOR PYRAMID %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tmp_disp(end-(sy2-sy1-1):end,:,:) = [];
    tmp_disp(:,end-(sx2-sx1-1):end,:) = [];
    D{frame}=tmp_disp;
end

% Preparing the output
Disparity = zeros(size(frames));
for frame = 1:n_frames
    % Horizontal disparity
    Disparity(:,:,frame,1)=D{frame}(:,:,1);
    % Vertical disparity
    Disparity(:,:,frame,2)=D{frame}(:,:,2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [II,x_pix,y_pix] = image_pyramid(II,n_scales)


[sy, sx, dummy] = size(II);

x_pix = cell(1,n_scales);
y_pix = cell(1,n_scales);

[ x_pix{n_scales}, y_pix{n_scales} ] = meshgrid(1:sx,1:sy);


lpf = [1 4 6 4 1]/16;

tmp = II;
II = cell(1,n_scales);
II{n_scales} = tmp;

for scale = n_scales-1:-1:1
    for frame = 1:2
        tmp(:,:,frame) = conv2b(conv2b(tmp(:,:,frame),lpf,3),lpf',3);
    end
    [Ny, Nx, dummy] = size(tmp);

    tmp = tmp(1:2:Ny,1:2:Nx,:);
    II{scale} = tmp;
    x_pix{scale} = x_pix{scale+1}(1:2:Ny,1:2:Nx);
    y_pix{scale} = y_pix{scale+1}(1:2:Ny,1:2:Nx);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = expand(D)
sy = size(D,1);
sx = size(D,2);

[X, Y] = meshgrid(1:(sx-1)/(2*sx-1):sx, ...
    1:(sy-1)/(2*sy-1):sy);

D = bilin_interp(D,X,Y);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [G]= warp_sequence(G,D)

[sy, sx, n_orient, phase_num, re_im] = size(G{2});

[X Y] = meshgrid(1:sx,1:sy);
clear sy sx n_orient phase_num re_im;

D(isnan(D)) = 0;

%% SLOW
G{2,1} = bilin_interp_orig(double(G{2,1}), X-D(:,:,1), Y-D(:,:,2));
G{2,2} = bilin_interp_orig(double(G{2,2}), X-D(:,:,1), Y-D(:,:,2));
%% FAST
% G{1,1} = bilin_interp_orig_mex(double(G{1,1}), X+D(:,:,1), Y+D(:,:,2));
% G{1,2} = bilin_interp_orig_mex(double(G{1,2}), X+D(:,:,1), Y+D(:,:,2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function D = merge_disparity(D1,D2)


invalid1 = isnan(D1);
invalid2 = isnan(D2);

D1(invalid1) = 0;
D2(invalid2) = 0;

invalid = invalid1 & invalid2;

% invalid = invalid1; % If unreliable at lower res, do not reattempt at
% higher res (reliability measure not suited for
% this). Low res valid estimates are retained even if
% they are never refined.

D = D1 + D2;
D(invalid) = NaN;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = population(G,energy_thr,ori_thr,ph_shift_file,scale)


[sy, sx, n_orient, phase_num, re_im] = size(G{2});

CL = G{1,1}; % REAL LEFT
CR = G{2,1}; % REAL RIGHT
SL = G{1,2}; % IMAG LEFT
SR = G{2,2}; % IMAG RIGHT

clear G;

%% ALLOCATE MEMORY FOR ENERGY
E=zeros(sy,sx,n_orient,phase_num);
Er=zeros(sy,sx,n_orient,phase_num);
El=zeros(sy,sx,n_orient);

%% ENERGY MODEL
for o = 1:n_orient
    for ph_n=1:phase_num
        
        % BINOCULAR ENERGY
        E(:,:,o,ph_n) = sqrt((CL(:,:,o) + squeeze(CR(:,:,o,ph_n))).^2 + (SL(:,:,o) + squeeze(SR(:,:,o,ph_n))).^2);
        E(:,:,o,ph_n)=squeeze(E(:,:,o,ph_n));  
        
        % MONOCULAR ENERGIES
        El(:,:,o)=(CL(:,:,o).^2 + squeeze(SL(:,:,o)).^2);
        Er(:,:,o,ph_n)=(SR(:,:,o,ph_n).^2 + squeeze(CR(:,:,o,ph_n)).^2);
        
    end
end

clear CL CR SL SR o ph_n;

%% THRESHOLDING 
% ENERGY THRESHOLD (computed on the monocular energy)

% ENERGY NORMALIZATION
El=repmat(El,[1 1 1 phase_num]);
El=El/max(max(max(max(El))));
Er=Er/max(max(max(max(Er))));

mask_energy=El>energy_thr & Er>energy_thr;

clear El Er;

% ORIENTATION THRESHOLD (VOID VALUES WHERE THE NUMBER OF VALID ORIENTATIONS IS UNDER THRESHOLD)   
mask_orientation=sum(sum(double(mask_energy),3)>=ori_thr,4)==phase_num;

% INTERSECTION OF ENERGY AND ORIENTATION THRESHOLD
mask=mask_energy.*repmat(mask_orientation,[1 1 n_orient phase_num]);

E=E.*mask;

clear mask mask_energy mask_orientation;

%% DECODING ALONG EACH ORIENTATION
E(isnan(E))=0;

% BINOCULAR ENERGY NORMALIZATION AND THRESHOLD (SCALE RESPONSE CURVES BETWEEN 0 AND 1)
E_max=repmat(max(E,[],4),[1 1 1 phase_num]);
E_min=repmat(min(E,[],4),[1 1 1 phase_num]);

E = (E-E_min)./(E_max-E_min);
clear E_max E_min;

E(E<0.7)=0;

% CENTER OF MASS DECODING ALONG THE ORIENTATIONS
load(ph_shift_file,'disp1');
Dori = sum(permute(repmat(disp1,[sy 1 sx n_orient]),[1 3 4 2]).*E,4)./sum(E,4);
clear disp1;

% COMPUTE X AND Y DISPARITY COMPONENTS FROM ORIENTATION COMPONENTS
theta=([1:n_orient]-1)*pi/8;
theta = permute(repmat(theta,[sy 1 sx]),[1 3 2]);

DX = Dori.*cos(theta);
DY = Dori.*sin(theta);

clear E theta;


%% INTERSECTION OF CONSTRAINT FOR VECTOR DISPARITY (GAUTAMA)
D = repmat(NaN, [sy sx 2]);
nc_min=3;

L_2 = DX.^2+DY.^2;
IND2 = isnan(L_2) | L_2<eps;

L_2(IND2) = Inf;
DX(IND2) = 0;
DY(IND2) = 0;
nc = sum(~IND2,3);

warning off

sumX = sum(DX,3);
sumY = sum(DY,3);
sumXYL_2 = sum(DX.*DY./L_2,3);
sumXXL_2 = sum(DX.^2./L_2,3);
sumYYL_2 = sum(DY.^2./L_2,3);
clear DX DY;
den = (sumXYL_2.^2-sumXXL_2.*sumYYL_2);

D(:,:,1)= -(sumX.*sumYYL_2-sumY.*sumXYL_2)./den;
D(:,:,2)= (sumX.*sumXYL_2-sumY.*sumXXL_2)./den;

warning on

D(repmat(nc,[1 1 2]) < nc_min) = NaN;


% %%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%
