%% Single filter plot
% Small script to visualise and test the filter reconstruction,
% notice that to speed up computation time, 2D filter convolution are 
% computed as 1D filters convolutions, therefore plotting the filters is 
% not straightforward but the filter needs to be reconstructed with the
% function gabor_filter.

addpath Functions

Orient=1; % Orientation
Phase=0; % Phase of the filter
Size='11x11'; % Dimension of the filter


filter=gabor_filter(Orient, Phase, Size);
subplot(1,2,1)
mesh(real(filter))
title('Even part of the filter')
subplot(1,2,2)
mesh(real(1i.*filter))
title('Odd part of the filter')

