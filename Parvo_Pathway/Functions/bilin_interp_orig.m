function I2 = bilin_interp_orig(I1,X,Y)

%
% function I2 = bilin_interp(I1,X,Y)
%
% Arbitrary image rewarping using bilinear interpolation
%
% X (column) and Y (row) (both floating point) are the source locations,
% used to fill the respective pixels
%
  
  
[nY1,nX1,rem] = size(I1); % source size
[nY2,nX2] = size(X); % target size

s1 = size(I1);
s = s1(3:end);

I2 = NaN.*zeros([ nY2 nX2 s ]);

for r = 1:rem

  for x = 1:nX2
    for y = 1:nY2
    
      % Pixel warping (2x2 group)
    
      x_w = floor(X(y,x));
      y_w = floor(Y(y,x));

      % Check validity

      if ( (x_w>0) & (x_w<nX1) & (y_w>0) & (y_w<nY1) )
        
        xs = X(y,x) - x_w;
        min_xs = 1-xs;
        ys = Y(y,x) - y_w;
        min_ys = 1-ys;
        
        w_00 = min_xs*min_ys;  % w_xy
        w_10 = xs*min_ys;
        w_01 = min_xs*ys;
        w_11 = xs*ys;
        
        I2(y,x,r) = w_00*I1(y_w,x_w,r) + w_10*I1(y_w,x_w+1,r) + ...
                  w_01*I1(y_w+1,x_w,r) + w_11*I1(y_w+1,x_w+1,r);

      end
    end
  end
  
end
