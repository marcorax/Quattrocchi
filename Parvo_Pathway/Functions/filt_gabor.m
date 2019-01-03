function [FI_n,f0] = filt_gabor(II,coarse_disp,fname,DC_thr)

if (nargin<3)
%   DC_thr = 1e-10;
DC_thr =0;
end


fltr = load(fname);
f0 = fltr.f0;
taps = fltr.taps;
F = fltr.F;
clear fltr;


n_orient = 8;

[nr,nc,n_frames,~] = size(II);
IF = zeros(nr,nc,n_orient,n_frames,2);

for camera = 1:2
    parfor frame = 1:n_frames
        tmp_result = zeros(nr,nc,n_orient);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Horizontal and vertical %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%


        % 0 orient

        F1Y = conv2b(II(:,:,frame,camera),F(1,:).',3);
        even = conv2b(F1Y,F(2,:),3);
        odd = conv2b(F1Y,F(3,:),3);

        tmp_result(:,:,1) = complex(even,odd);

        % pi/2 orient

        F1X = conv2b(II(:,:,frame,camera),F(1,:),3);
        even = conv2b(F1X,F(2,:).',3);
        odd = conv2b(F1X,F(3,:).',3);

        tmp_result(:,:,5) = complex(even,odd);


        %%%%%%%%%%%%
        % Diagonal %
        %%%%%%%%%%%%


        F4Y = conv2b(II(:,:,frame,camera),F(4,:).',3);

        F4YF4X = conv2b(F4Y,F(4,:),3);
        F4YF5X = conv2b(F4Y,F(5,:),3);

        F5Y = conv2b(II(:,:,frame,camera),F(5,:).',3);

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


        F8Y = conv2b(II(:,:,frame,camera),F(8,:).',3);

        F8YF6X = conv2b(F8Y,F(6,:),3);
        F8YF7X = conv2b(F8Y,F(7,:),3);

        F9Y = conv2b(II(:,:,frame,camera),F(9,:).',3);

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

        F6Y = conv2b(II(:,:,frame,camera),F(6,:).',3);

        F6YF8X = conv2b(F6Y,F(8,:),3);
        F6YF9X = conv2b(F6Y,F(9,:),3);

        F7Y = conv2b(II(:,:,frame,camera),F(7,:).',3);

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

        IF(:,:,:,frame,camera)=tmp_result;
    end
end

invalid = (abs(real(IF))<DC_thr) | (abs(imag(IF))<DC_thr);
IF(invalid) = NaN;

    FI_n{1}=real(IF);
    FI_n{2}=imag(IF);

