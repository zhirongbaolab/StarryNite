%Fernando Amat September 17th 2010
%Modification to allow sigma having different values

%If X is of type double, the mex file is faster. Otherwise this m-file is
%faster for single and int images
function I=imgaussianAnisotropy(I,sigma,siz)
% IMGAUSSIAN filters an 1D, 2D color/greyscale or 3D image with an 
% Gaussian filter. This function uses for filtering IMFILTER or if 
% compiled the fast  mex code imgaussian.c . Instead of using a 
% multidimensional gaussian kernel, it uses the fact that a Gaussian 
% filter can be separated in 1D gaussian kernels.
%
% J=IMGAUSSIAN(I,SIGMA,SIZE)
%
% inputs,
%   I: The 1D, 2D greyscale/color, or 3D input image with 
%           data type Single or Double
%   SIGMA: The sigma used for the Gaussian kernel
%   SIZE: Kernel size (single value) (default: sigma*6)
% 
% outputs,
%   J: The gaussian filtered image
%
% note, compile the code with: mex imgaussian.c -v
%
% example,
%   I = im2double(imread('peppers.png'));
%   figure, imshow(imgaussian(I,10));
% 
% Function is written by D.Kroon University of Twente (September 2009)

if(~exist('siz','var')), siz=sigma*6; end

if(length(sigma)~=length(size(I)))
    error 'You must specify one sigma for each dimension of the image'
end


% Filter each dimension with the 1D Gaussian kernels\
if(ndims(I)==1)
    % Make 1D Gaussian kernel
    kk=1;
    x=-ceil(siz(kk)/2):ceil(siz(kk)/2);
    H = exp(-(x.^2/(2*sigma(kk)^2)));
    H = H/sum(H(:));

    I=imfilter(I,H, 'same' ,'replicate');
elseif(ndims(I)==2)
    % Make 1D Gaussian kernel
    kk=1;
    x=-ceil(siz(kk)/2):ceil(siz(kk)/2);
    H = exp(-(x.^2/(2*sigma(kk)^2)));
    H = H/sum(H(:));
    Hx=reshape(H,[length(H) 1]);
    
    
    kk=2;
    x=-ceil(siz(kk)/2):ceil(siz(kk)/2);
    H = exp(-(x.^2/(2*sigma(kk)^2)));
    H = H/sum(H(:));
    Hy=reshape(H,[1 length(H)]);
    I=imfilter(imfilter(I,Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate');
elseif(ndims(I)==3)
    
    
    if(size(I,3)<4) % Detect if 3D or color image
        % Make 1D Gaussian kernel
        kk=1;
        x=-ceil(siz(kk)/2):ceil(siz(kk)/2);
        H = exp(-(x.^2/(2*sigma(kk)^2)));
        H = H/sum(H(:));
        Hx=reshape(H,[length(H) 1]);
        
        % Make 1D Gaussian kernel
        kk=2;
        x=-ceil(siz(kk)/2):ceil(siz(kk)/2);
        H = exp(-(x.^2/(2*sigma(kk)^2)));
        H = H/sum(H(:));
        Hy=reshape(H,[1 length(H)]);
        for k=1:size(I,3)
            I(:,:,k)=imfilter(imfilter(I(:,:,k),Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate');
        end
    else
        % Make 1D Gaussian kernel
        kk=1;
        x=-ceil(siz(kk)/2):ceil(siz(kk)/2);
        H = exp(-(x.^2/(2*sigma(kk)^2)));
        H = H/sum(H(:));
        Hx=reshape(H,[length(H) 1 1]);
        
        % Make 1D Gaussian kernel
        kk=2;
        x=-ceil(siz(kk)/2):ceil(siz(kk)/2);
        H = exp(-(x.^2/(2*sigma(kk)^2)));
        H = H/sum(H(:));
        Hy=reshape(H,[1 length(H) 1]);
        
        % Make 1D Gaussian kernel
        kk=3;
        x=-ceil(siz(kk)/2):ceil(siz(kk)/2);
        H = exp(-(x.^2/(2*sigma(kk)^2)));
        H = H/sum(H(:));
        Hz=reshape(H,[1 1 length(H)]);
        
        I=imfilter(imfilter(imfilter(I,Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate'),Hz, 'same' ,'replicate');
    end
else
    error('imgaussian:input','unsupported input dimension');
end
