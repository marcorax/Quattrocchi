function g=gabor_filter(Orient, Phase, Size)
%% gabor_filter(Orient,Phase,Size)
% Function that allows to create a simple filter for visualisation purpouses
% Orient: is an integer value selecting the orientation of the filter,
%         0=0 1=pi/2 2=pi/4 3=3pi/4 4=pi/8 5=7pi/8 6=3pi/8 7=5/pi8
% Phase: is the phase of the filter 
% Size: is a sting selecting the dimension of the filters, 
%       only two values are accepted '11x11' '34x34'

if(strcmp(Size,'11x11'))
    fname='Filters/Gt11B0.0833f0.25.mat';
    ph_shift_file='Filters/ph_shift_PH7_F0.25.mat';
end

if(strcmp(Size,'34x34'))
    fname='Filters/Gt43B0.0208f0.063.mat';
    ph_shift_file='Filters/ph_shift_PH7_F0.063.mat';
end

fltr = load(fname);
F = fltr.F;
clear fltr;
switch Orient
    case 0
        f=complex((F(1,:).'*F(2,:)),(F(1,:).'*F(3,:)));
    case 1
        f=complex((F(2,:).'*F(1,:)),(F(3,:).'*F(1,:)));
    case 2
        f=complex((F(4,:).'*F(4,:))-(F(5,:).'*F(5,:)),(F(5,:).'*F(4,:))+(F(4,:).'*F(5,:)));
    case 3
        f=complex((F(4,:).'*F(4,:))+(F(5,:).'*F(5,:)),(F(5,:).'*F(4,:))-(F(4,:).'*F(5,:))); 
    case 4
        f=complex((F(8,:).'*F(6,:))-(F(9,:).'*F(7,:)),(F(9,:).'*F(6,:))+(F(8,:).'*F(7,:)));
    case 5
        f=complex((F(8,:).'*F(6,:))+(F(9,:).'*F(7,:)),(F(9,:).'*F(6,:))-(F(8,:).'*F(7,:)));
    case 6
        f=complex((F(6,:).'*F(8,:))-(F(7,:).'*F(9,:)),(F(7,:).'*F(8,:))+(F(6,:).'*F(9,:)));
    case 7
        f=complex((F(6,:).'*F(8,:))+(F(7,:).'*F(9,:)),(F(7,:).'*F(8,:))-(F(6,:).'*F(9,:)));
end
        
        
if Phase~=0
load(ph_shift_file) 
G_even_tmp=real(f);
G_odd_tmp=real(1i*f);

ph=ph_shift(1,Phase);
G_even=G_even_tmp*cos(ph)-G_odd_tmp*sin(ph);
G_odd=G_odd_tmp*cos(ph)+G_even_tmp*sin(ph);
g=G_even+1i.*G_odd;

else
    g=f;
    
end
