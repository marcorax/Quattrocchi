function IF = filt_sep_gabor(I,fname,DC_thr)

% function IF = filt_sep_gabor(I,fname,DC_thr)
%  

if (nargin<3)
  DC_thr = 1e-10;
end


fltr = load(fname);
taps = fltr.taps;
F = fltr.F;
clear fltr;


n_orient = 8;

[nr,nc,n_frames] = size(I);
padsize = floor(taps/2); 
IF = zeros(nr,nc,n_orient,n_frames);

for frame = 1:n_frames


  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Horizontal and vertical %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  % 0

  F1Y = conv2b(I(:,:,frame),F(1,:).',3);
  even = conv2b(F1Y,F(2,:),3);
  odd = conv2b(F1Y,F(3,:),3);
  
  IF(:,:,1,frame) = complex(even,odd);

  clear even odd F1Y;
  
    % pi

  F1Y = conv2b(I(:,:,frame),F(1,:).',3);
  even = conv2b(F1Y,F(2,:),3);
  odd = -(conv2b(F1Y,F(3,:),3));
  
  IF(:,:,9,frame) = complex(even,odd);

  clear even odd F1Y;
  
  % pi/2

  F1X = conv2b(I(:,:,frame),F(1,:),3);
  even = conv2b(F1X,F(2,:).',3);
  odd = conv2b(F1X,F(3,:).',3);

  IF(:,:,5,frame) = complex(even,odd);

  clear even odd F1X;
  
  % 11 pi/8
  
  F1X = conv2b(I(:,:,frame),F(1,:),3);
  even = conv2b(F1X,F(2,:).',3);
  odd = -conv2b(F1X,F(3,:).',3);

  IF(:,:,13,frame) = complex(even,odd);

  clear even odd F1X;

  
  %%%%%%%%%%%%
  % Diagonal %
  %%%%%%%%%%%%
  
  
  F4Y = conv2b(I(:,:,frame),F(4,:).',3);

  F4YF4X = conv2b(F4Y,F(4,:),3);
  F4YF5X = conv2b(F4Y,F(5,:),3);

  clear F4Y;
  
  F5Y = conv2b(I(:,:,frame),F(5,:).',3);

  F5YF4X = conv2b(F5Y,F(4,:),3);
  F5YF5X = conv2b(F5Y,F(5,:),3);

  clear F5Y;
  
  % pi/4

  even = F4YF4X - F5YF5X;
  odd = F5YF4X + F4YF5X;
  
  IF(:,:,3,frame) = complex(even,odd);
  
    % 10 pi/8

  even = F4YF4X - F5YF5X;
  odd = -(F5YF4X + F4YF5X);
  
  IF(:,:,11,frame) = complex(even,odd);
  
  % 3pi/4

  even = F4YF4X + F5YF5X;
  odd = F5YF4X - F4YF5X;
  
  IF(:,:,7,frame) = complex(even,odd);


  
    % 14pi/4

  even = F4YF4X + F5YF5X;
  odd = -(F5YF4X - F4YF5X);
  
  IF(:,:,15,frame) = complex(even,odd);

  clear even odd F4YF4X F4YF5X F5YF4X F5YF5X;

  
  %%%%%%%%%%%%%%%%%
  % 'In-betweens' %
  %%%%%%%%%%%%%%%%%

  
  F8Y = conv2b(I(:,:,frame),F(8,:).',3);
  
  F8YF6X = conv2b(F8Y,F(6,:),3);
  F8YF7X = conv2b(F8Y,F(7,:),3);

  clear F8Y;
  
  F9Y = conv2b(I(:,:,frame),F(9,:).',3);

  F9YF6X = conv2b(F9Y,F(6,:),3);
  F9YF7X = conv2b(F9Y,F(7,:),3);
  
  clear F9Y;

  % pi/8
  
  even = F8YF6X - F9YF7X;
  odd = F9YF6X + F8YF7X;
  
  IF(:,:,2,frame) = complex(even,odd);
  
    % 9pi/8
  
  even = F8YF6X - F9YF7X;
  odd = -(F9YF6X + F8YF7X);
  
  IF(:,:,10,frame) = complex(even,odd);
  
  % 7pi/8
  
  even = F8YF6X + F9YF7X;
  odd = F9YF6X - F8YF7X;
  
  IF(:,:,8,frame) = complex(even,odd);
  
    % 15pi/8
  
  even = F8YF6X + F9YF7X;
  odd = -(F9YF6X - F8YF7X);
  
  IF(:,:,16,frame) = complex(even,odd);

  clear even odd F8YF6X F8YF7X F9YF6X F9YF7X;
  
  
  F6Y = conv2b(I(:,:,frame),F(6,:).',3);
  
  F6YF8X = conv2b(F6Y,F(8,:),3);
  F6YF9X = conv2b(F6Y,F(9,:),3);

  clear F6Y;
  
  F7Y = conv2b(I(:,:,frame),F(7,:).',3);

  F7YF8X = conv2b(F7Y,F(8,:),3);
  F7YF9X = conv2b(F7Y,F(9,:),3);
  
  clear F7Y;

  % 3pi/8
  
  even = F6YF8X - F7YF9X;
  odd = F7YF8X + F6YF9X;
  
  IF(:,:,4,frame) = complex(even,odd);
  
    % 11pi/8
  
  even = F6YF8X - F7YF9X;
  odd = -(F7YF8X + F6YF9X);
  
  IF(:,:,12,frame) = complex(even,odd);

  % 5pi/8
  
  even = F6YF8X + F7YF9X;
  odd = F7YF8X - F6YF9X;
  
  IF(:,:,6,frame) = complex(even,odd);
  
    % 13pi/8
  
  even = F6YF8X + F7YF9X;
  odd = -(F7YF8X - F6YF9X);
  
  IF(:,:,14,frame) = complex(even,odd);
  
  clear even odd F6YF8X F6YF9X F7YF8X F7YF9X;
 
end



invalid = (abs(real(IF))<DC_thr) | (abs(imag(IF))<DC_thr);
IF(invalid) = NaN;

G_even=real(IF);
G_odd=imag(IF);
clear IF;
IF{1}=G_even;
IF{2}=G_odd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%