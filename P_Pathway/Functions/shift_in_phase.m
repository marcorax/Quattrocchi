

function G=shift_in_phase(F_real,F_imag,file_name)

n_filters=size(F_real,3);
load(file_name)
[n_orient n_phase]=size(ph_shift);
%switch n_filters
%     case 8
    p=ph_shift;
%     case 8
%     p=ph_shift(1:2:end,:);
%     n_orient=n_filters;
%     case 4
%     p=ph_shift(1:2:7,:);
%     n_orient=n_filters;
%     case 2
%     p=ph_shift(3:4:7,:);
%     n_orient=n_filters;
%end
for orient = 1:n_orient
    G_even_tmp=F_real(:,:,orient);
    G_odd_tmp=F_imag(:,:,orient);
    phase=1;

    for(ph=p(orient,:))
        G_even(:,:,orient,phase)=G_even_tmp*cos(ph)-G_odd_tmp*sin(ph);
        G_odd(:,:,orient,phase)=G_odd_tmp*cos(ph)+G_even_tmp*sin(ph);
        phase=phase+1;
    end
end
G{1,1} = G_even;
G{1,2} = G_odd;

