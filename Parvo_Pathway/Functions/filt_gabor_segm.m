function [FI_n,f0] = filt_gabor_segm(II,seg_masks,coarse_disparity,fname,DC_thr)

if (nargin<5)
%   DC_thr = 1e-10;
DC_thr =0;
end


fltr = load(fname);
f0 = fltr.f0;
F = fltr.F;
clear fltr;

n_orient = 8;

[nr,nc,n_frames,~] = size(II);
IF = zeros(nr,nc,n_orient,n_frames,2);

for camera = 1:2
    parfor frame = 1:n_frames
        
        tmp_result = zeros(nr,nc,n_orient);
        current_mask = seg_masks(:,:,frame,camera);
        current_frame = II(:,:,frame,camera).*current_mask;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Horizontal and vertical %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%                
        
        % 0 orient

        F1Y = conv2b(current_frame,F(1,:).',3);
        even = conv2b(F1Y,F(2,:),3);
        odd = conv2b(F1Y,F(3,:),3);

        tmp_result(:,:,1) = complex(even,odd);

        % pi/2 orient

        F1X = conv2b(current_frame,F(1,:),3);
        even = conv2b(F1X,F(2,:).',3);
        odd = conv2b(F1X,F(3,:).',3);

        tmp_result(:,:,5) = complex(even,odd);


        %%%%%%%%%%%%
        % Diagonal %
        %%%%%%%%%%%%


        F4Y = conv2b(current_frame,F(4,:).',3);

        F4YF4X = conv2b(F4Y,F(4,:),3);
        F4YF5X = conv2b(F4Y,F(5,:),3);

        F5Y = conv2b(current_frame,F(5,:).',3);

        F5YF4X = conv2b(F5Y,F(4,:),3);
        F5YF5X = conv2b(F5Y,F(5,:),3);

        % pi/4 orient

        even = F4YF4X - F5YF5X;
        odd = F5YF4X + F4YF5X;

        tmp_result(:,:,3) = complex(even,odd);

        % 3pi/4 orient

        even = F4YF4X + F5YF5X;
        odd = F5YF4X - F4YF5X;

        tmp_result(:,:,7) = complex(even,odd);


        %%%%%%%%%%%%%%%%%
        % 'In-betweens' %
        %%%%%%%%%%%%%%%%%


        F8Y = conv2b(current_frame,F(8,:).',3);

        F8YF6X = conv2b(F8Y,F(6,:),3);
        F8YF7X = conv2b(F8Y,F(7,:),3);

        F9Y = conv2b(current_frame,F(9,:).',3);

        F9YF6X = conv2b(F9Y,F(6,:),3);
        F9YF7X = conv2b(F9Y,F(7,:),3);

        % pi/8 orient

        even = F8YF6X - F9YF7X;
        odd = F9YF6X + F8YF7X;

        tmp_result(:,:,2) = complex(even,odd);

        % 7pi/8 orient

        even = F8YF6X + F9YF7X;
        odd = F9YF6X - F8YF7X;

        tmp_result(:,:,8) = complex(even,odd);  

        F6Y = conv2b(current_frame,F(6,:).',3);

        F6YF8X = conv2b(F6Y,F(8,:),3);
        F6YF9X = conv2b(F6Y,F(9,:),3);

        F7Y = conv2b(current_frame,F(7,:).',3);

        F7YF8X = conv2b(F7Y,F(8,:),3);
        F7YF9X = conv2b(F7Y,F(9,:),3);

        % 3pi/8 orient

        even = F6YF8X - F7YF9X;
        odd = F7YF8X + F6YF9X;

        tmp_result(:,:,4) = complex(even,odd);

        % 5pi/8

        even = F6YF8X + F7YF9X;
        odd = F7YF8X - F6YF9X;

        tmp_result(:,:,6) = complex(even,odd);
        tmp_result = position_shift(tmp_result,squeeze(coarse_disparity(:,:,frame,:)));
        IF(:,:,:,frame,camera)=tmp_result;
    end
end

invalid = (abs(real(IF))<DC_thr) | (abs(imag(IF))<DC_thr);
IF(invalid) = NaN;

    FI_n{1}=real(IF);
    FI_n{2}=imag(IF);
    
function result = position_shift(matrix, disparity)
% This function shifts the results of gabor filtering with previously
% computed disparity, in a way that each parvo disparity neuron ad the end
% will inherit position shift from the magno disparity neurons
% matrix is a 3D matrix where the first two dimensions are the same as 
% the cut frames and the third dimension is the orientation of the filter
% Disaparity is a 3D matrix where the two first dimension are the same as 
% the cut frames and the third one is 0 for horizonal disparity and 1 for
% vertical disparity
[ly,lx,~]=size(matrix);
result = matrix;
[indy,indx] = find(abs(disparity(:,:,1))+abs(disparity(:,:,2))~=0);
indy_shift = int32(indy+round(disparity(sub2ind(size(disparity),indy,indx,ones(size(indy))*2))));
indx_shift = int32(indx+round(disparity(sub2ind(size(disparity),indy,indx,ones(size(indy))*1))));
% remove all shifts outside boundaries
indy_shift = indy_shift((0>indy_shift)&(indy_shift>=ly));
indx_shift = indx_shift((0>indx_shift)&(indx_shift>=lx));
indy = indy((0>indy_shift)&(indy_shift>=ly));
indx = indx((0>indx_shift)&(indx_shift>=lx));
for orient = 1:8
result(sub2ind(size(matrix),indy,indx,ones(size(indy))*orient)) = ...
    matrix(sub2ind(size(matrix),indy_shift,indx_shift,ones(size(indy))*orient));
end

