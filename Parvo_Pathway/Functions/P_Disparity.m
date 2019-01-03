function D = P_Disparity(II,coarse_disp,energy_thr,ori_thr,ph_shift_file,filter_file)

% TODO, use the convolution on the segmented part of the image
% Compute disparity refining with gabor filter 
[FI,~] = filt_gabor(II,coarse_disp,filter_file);

% TODO modify the code to accept position shift 
% Compute phase shift
G{1,1}=FI{1}(:,:,:,1); G{1,2}=FI{2}(:,:,:,1);
G(2,:)=shift_in_phase(FI{1}(:,:,:,2),FI{2}(:,:,:,2),ph_shift_file);
clear FI

% Plot Filter shape (uncomment to use it)
% plot_filter(G)

% Population decoding
D = population(G,energy_thr,ori_thr,ph_shift_file);

% Mask unreliable values
invalid = isnan(D);
D(invalid) = NaN;

clear E G
 


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
theta=(1:n_orient-1)*pi/8;
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