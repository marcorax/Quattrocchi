function imf=conv2b(im, ker, border);
% CONV2B Two-dimensional convolution with border handling
%	Performs 2-D convolution returning an image with the same 
%	dimensions as the input image.
%	The border handling is specified by border:
%		0 -> pad with zeros
%		1 -> periodic convolution
%		2 -> specular borders
%       3 -> repeat edge borders
%
%	imf=conv2b(im, ker, border);
%
% This function is not yet prepared for even sized kernels
%
% ONG - 10/96	
%


% checking dimensions
[nky, nkx]=size(ker);
if (~rem(nkx,2))||(~rem(nky,2))
	fprintf(1,'\n\t CONV2B WARNING: not specific handling for even sized kernels\n\n');
end

% proper treatment for each value of border
if border==0
	sh='same';
elseif (border==1)||(border==2)||(border==3)
	sh='valid';
	Bx=(nkx-1)/2;
	By=(nky-1)/2;
	im=putborde(im,Bx,By,border);
else
	error('Not a valid value for border')
end

% convolving the image
imf=conv2(im,ker,sh);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function imb=putborde(im,Nx,Ny,method);
% PUTBORDE Puts a border to an image:
%	method=1 -> consider the image periodical
%	method=2 -> specular
%	method=3 -> repeat the edge pixel
%	
%	imb=putborde(im,Nx,Ny,method);
%
% ONG - 10/96
%

[sy, sx]=size(im);
imb=zeros(sy+2*Ny,sx+2*Nx);
imb(1+Ny:sy+Ny,1+Nx:sx+Nx)=im;

if method == 1
	imb(1:Ny,1+Nx:sx+Nx)=im(sy-Ny+1:sy,:);
	imb(sy+Ny+1:sy+2*Ny,1+Nx:sx+Nx)=im(1:Ny,:);
	imb(1+Ny:sy+Ny,1:Nx)=im(:,sx-Nx+1:sx);
	imb(1+Ny:sy+Ny,sx+Nx+1:sx+2*Nx)=im(:,1:Nx);
	imb(1:Ny,1:Nx)=im(sy-Ny+1:sy,sx-Nx+1:sx);
	imb(1:Ny,sx+Nx+1:sx+2*Nx)=im(sy-Ny+1:sy,1:Nx);
	imb(sy+Ny+1:sy+2*Ny,1:Nx)=im(1:Ny,sx-Nx+1:sx);
	imb(sy+Ny+1:sy+2*Ny,sx+Nx+1:sx+2*Nx)=im(1:Ny,1:Nx);
elseif method == 2
	imb(Ny:-1:1,1+Nx:sx+Nx)=im(2:Ny+1,:);
	imb(sy+2*Ny:-1:sy+Ny+1,1+Nx:sx+Nx)=im(sy-Ny:sy-1,:);
	imb(1+Ny:sy+Ny,Nx:-1:1)=im(:,2:Nx+1);
	imb(1+Ny:sy+Ny,sx+2*Nx:-1:sx+Nx+1)=im(:,sx-Nx:sx-1);
	imb(Ny:-1:1,Nx:-1:1)=im(2:Ny+1,2:Nx+1);
	imb(Ny:-1:1,sx+2*Nx:-1:sx+Nx+1)=im(2:Ny+1,sx-Nx:sx-1);
	imb(sy+2*Ny:-1:sy+Ny+1,Nx:-1:1)=im(sy-Ny:sy-1,2:Nx+1);
	imb(sy+2*Ny:-1:sy+Ny+1,sx+2*Nx:-1:sx+Nx+1)=im(sy-Ny:sy-1,sx-Nx:sx-1);
elseif method==3
	for k=1:Nx
		imb(Ny+1:sy+Ny,k)=im(:,1);
		imb(Ny+1:sy+Ny,k+sx+Nx)=im(:,sx);
	end
	for k=1:Ny
		imb(k,Nx+1:sx+Nx)=im(1,:);
		imb(k+sy+Ny,Nx+1:sx+Nx)=im(sy,:);
	end
else
	error('Not a valid value for method')
end

