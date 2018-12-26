function plot_filter(G)


RealRight = G{2,1};
ImagRight = G{2,2};

n_orient = size(RealRight,3);
phase_num = size(RealRight,4);

figure

% for j=0:n_orient-1
j=0;
    for i=0:phase_num-1
        
        subplot(3,3,(i+1)+j*phase_num)
        imagesc(RealRight(:,:,j+1,i+1)), axis square, colormap(gray)
    end
% end
text(-120,-160,'Filter real part')
text(-100,25,'Phase')
text(-220,-60,'Orient')

figure

% for j=0:n_orient-1
    for i=0:phase_num-1
        
        subplot(3,3,(i+1)+j*phase_num)
        imagesc(ImagRight(:,:,j+1,i+1)), axis square, colormap(gray)
    end
% end
text(-120,-160,'Filter imaginary part')
text(-100,25,'Phase')
text(-220,-60,'Orient')
