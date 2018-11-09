function [F,f0] = filt_steer_or_gabor(I,fname,DC_thr)

% function [F,f0] = filt_steer_or_gabor(I,fname,DC_thr)


if (nargin<3)
%   DC_thr = 1e-10;
DC_thr =0;
end


fltr = load(fname);

switch fltr.type
  
 case 's'

  F = filt_sep_steer(I,fname,DC_thr);
  
 case 'g'

  F = filt_sep_gabor(I,fname,DC_thr);
  
end

f0 = fltr.f0;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function F = filt_sep_steer(I,fname,DC_thr)

% function F = filt_sep_steer(I,fname,DC_thr)
%  

fltr = load(fname);
taps = fltr.taps;
G = fltr.G;
H = fltr.H;
order = fltr.order;
clear fltr;

n_orient = 8;

[nr,nc,n_frames] = size(I);
  
nG = size(G,1);
nH = size(H,1);

F = zeros(nr,nc,n_orient,n_frames);

for frame = 1:n_frames

  % Construct the basis filter responses

  Gb = zeros(nr,nc,nG);
  Hb = zeros(nr,nc,nH);

  for i = 1:nG

    % X-direction
  
    f = squeeze(G(i,1,:))';
    T = conv2b(I(:,:,frame),f,3);

    % Y-direction
  
    f = squeeze(G(i,2,:));
    Gb(:,:,i) = conv2b(T,f,3);
  
  end

  for i = 1:nH

    % X-direction
  
    f = squeeze(H(i,1,:))';
    T = conv2b(I(:,:,frame),f,3);

    % Y-direction
  
    f = squeeze(H(i,2,:));
    Hb(:,:,i) = conv2b(T,f,3);
    
  end
  
  % Steer filters and compute oriented responses

  padsize = floor(taps/2); 
  
  angles = -(0:pi/n_orient:pi-(pi/n_orient)); % negate because gabors
                                              % rotate in the other direction 

  for o = 1:n_orient

    even = zeros(nr,nc);
    odd = zeros(nr,nc);
    
    [Gk,Hk] = gen_steer_coeff(angles(o),order);
    
    for i = 1:nG
      even = even + Gk(i).*Gb(:,:,i);
    end
    
    for i = 1:nH
      odd = odd + Hk(i).*Hb(:,:,i);
    end
    
    invalid = (abs(even)<DC_thr) | (abs(odd)<DC_thr);
    
    even(invalid) = NaN;
    odd(invalid) = NaN;
    
    F(:,:,o,frame) = complex(even,odd);
    
  end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function FI_n = filt_sep_gabor(I,fname,DC_thr)

% function IF = filt_sep_gabor(I,fname,DC_thr)
%  

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
  
  % pi/2

  F1X = conv2b(I(:,:,frame),F(1,:),3);
  even = conv2b(F1X,F(2,:).',3);
  odd = conv2b(F1X,F(3,:).',3);

  IF(:,:,5,frame) = complex(even,odd);

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
  
  % 3pi/4

  even = F4YF4X + F5YF5X;
  odd = F5YF4X - F4YF5X;
  
  IF(:,:,7,frame) = complex(even,odd);

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
  
  % 7pi/8
  
  even = F8YF6X + F9YF7X;
  odd = F9YF6X - F8YF7X;
  
  IF(:,:,8,frame) = complex(even,odd);

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

  % 5pi/8
  
  even = F6YF8X + F7YF9X;
  odd = F7YF8X - F6YF9X;
  
  IF(:,:,6,frame) = complex(even,odd);
  
  clear even odd F6YF8X F6YF9X F7YF8X F7YF9X;
 
end



invalid = (abs(real(IF))<DC_thr) | (abs(imag(IF))<DC_thr);
IF(invalid) = NaN;

    FI_n{1}=real(IF);
    FI_n{2}=imag(IF);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Gk,Hk] = gen_steer_coeff(angle,order)

% function [Gk,Hk] = gen_steer_coeff(angle,order)
% 
% angle in radians
  

if (nargin<2)
  order = 4;
end

if ~ismember(order,[2 4])
  error('Order should be equal to 2 or 4.\n');
end

c = cos(angle);
s = sin(angle);

switch order

  
 case 2 % Second order filters

  
  Gk = zeros(3,1);
  Hk = zeros(4,1);

  Gk = [ c.^2
         -2.*c.*s
         s.^2 ];
  
  Hk = [ c.^3
         -3.*c.^2.*s
         3.*c.*s.^2
         -s.^3 ];

  
 case 4 % Fourth order filters

  
  Gk = zeros(5,1);
  Hk = zeros(6,1);

  Gk = [ c.^4
         -4.*c.^3.*s
         6.*c.^2.*s.^2
         -4.*c.*s.^3
         s.^4 ];

  Hk = [ c.^5
         -5.*c.^4.*s
         10.*c.^3.*s.^2
         -10.*c.^2.*s.^3
         5.*c.*s.^4
         -s.^5 ];
  
end
