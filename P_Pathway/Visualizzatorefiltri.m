%% Visualizzatore di filtri orientati

% 0 = 0 
% 1 = pi/2 
% 2 = pi/4 
% 3 = 3pi/4 
% 4 = pi/8 
% 5 = 7pi/8 
% 6 = 3pi/8 
% 7 = 5/pi8

i=1;
k=0;

f=gaborfilter(i,k);
subplot(1,2,1)
mesh(real(f))
title('Even')
subplot(1,2,2)
mesh(real(1i.*f))
title('Odd')

