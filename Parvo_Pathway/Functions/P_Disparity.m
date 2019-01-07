function fine_disp = P_Disparity(II,coarse_disp,seg_masks,energy_thr,ori_thr,ph_shift_file,filter_file)

% Put NaN seg_masks to 0 for computing the convolution
seg_masks(isnan(seg_masks))=0;

% Compute disparity refining with gabor filter 
[FI,~] = filt_gabor_segm(II,seg_masks,coarse_disp,filter_file);

[~,~,n_frames,~]=size(II);
D = cell(length(n_frames),1);
parfor frame = 1:n_frames 
    G = cell(2);
    current_FI_real = squeeze(FI{1}(:,:,:,frame,:));
    current_FI_imag = squeeze(FI{2}(:,:,:,frame,:));

    % Compute phase shift
    G{1,1}=current_FI_real(:,:,:,1); G{1,2}=current_FI_imag(:,:,:,1);
    G(2,:)=shift_in_phase(current_FI_real(:,:,:,2),current_FI_imag(:,:,:,2),ph_shift_file);

    % Population decoding
    D{frame} = population(G,energy_thr,ori_thr,ph_shift_file);

    % Mask unreliable values
    invalid = isnan(D{frame});
    D{frame}(invalid) = NaN;
end
% Preparing the output
fine_disp = zeros(size(coarse_disp));
for frame = 1:n_frames
    % Horizontal disparity
    fine_disp(:,:,frame,1)=D{frame}(:,:,1);
    % Vertical disparity
    fine_disp(:,:,frame,2)=D{frame}(:,:,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = population(G,energy_thr,ori_thr,ph_shift_file)


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
theta=(0:n_orient-1)*pi/8;
theta = permute(repmat(theta,[sy 1 sx]),[1 3 2]);

DX = Dori.*cos(theta);
DY = Dori.*sin(theta);

clear E theta;


%% INTERSECTION OF CONSTRAINT FOR VECTOR DISPARITY (GAUTAMA)
D = NaN([sy sx 2]);
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