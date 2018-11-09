function g=gaborfilter(i,theta,fname,ph_shift_file)
%% gaborfilter(i,k,fname)
% Funzione che ottiene un filtro dal set 9*11 del prof
% i permette di scegliere l'orientamento del filtro,
% 0=0 1=pi/2 2=pi/4 3=3pi/4 4=pi/8 5=7pi/8 6=3pi/8 7=5/pi8
% theta è l'indice che indica l'orientamento estratto 
% da ph_shift_file caricato
% fname è la stringa della posizione dei filtri
% se non specificata il path di default
% è 'FILTERS/Gt11B0.0833f0.25.mat'
% Così come per ph_shift_file la cui path di default è
% 'FILTERS\ph_shift_PH7_F0.25.mat'

if(nargin<3)
    fname='FILTERS/Gt11B0.0833f0.25.mat';
end

if(nargin<4)
    ph_shift_file='FILTERS\ph_shift_PH7_F0.25.mat';
end

fltr = load(fname);
F = fltr.F;
clear fltr;
switch i
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
        
        
if theta~=0
load(ph_shift_file)
G_even_tmp=real(f);
G_odd_tmp=real(1i*f);

ph=ph_shift(1,theta);
G_even=G_even_tmp*cos(ph)-G_odd_tmp*sin(ph);
G_odd=G_odd_tmp*cos(ph)+G_even_tmp*sin(ph);
g=G_even+1i.*G_odd;

else
    g=f;
    
end
