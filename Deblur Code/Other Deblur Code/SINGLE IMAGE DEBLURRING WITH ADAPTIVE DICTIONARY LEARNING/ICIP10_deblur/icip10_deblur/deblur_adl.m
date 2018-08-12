function [deblur_im,deblur_kernel] = deblur_adl(ori_im,kernel_size,patch_size,iter_num,loc)
% function deblur_adl is based on the paper 
% <Single image deblurring with adaptive dictionary learning> 
% It tries to solve the single image deblur problem using sparse prior and 
% adaptive dictionary.
% 
% Required fields in PARAMS:
%  --------------------------
%    'ori_im' - blurry image
% 
% Optional fields in PARAMS:
%  --------------------------
%    'kernel_size' - size of blur kernel
%    'patch_size' - size of patch used for dictionary learning and sparse 
%                  representation, related to dictionary basis size
%    'iter_num' - to control the number of iterations
%    'loc' - specify the image window used to estimate kernel, if the
%           window size is large, the speed will be slow, better no more than
%           200x200
%           
%
% Zhe Hu
% Computer Science Department
% University of California, Merced
% zhu@ucmerced.edu
% 
% October 2009

[row,col,channel] = size(ori_im);

if nargin<2
    kernel_size = 7;
end
if nargin<3
    patch_size = 12;
end
if nargin<4
   iter_num = 6;
end
if nargin<5
   loc =  [1,row,1,col];
end

% choose DCT as initial dictionary
dictionary_size = 10*patch_size^2;
initdict = odctdict(patch_size^2, dictionary_size);

%% read image(or image patch) 
Blurred = im2double(ori_im);
if isrgb(ori_im)
   ori_im = rgb2gray(ori_im);
end
Blurred_gray = im2double(Blurred);

%% extract the subwindow for deblurring
Blurred_win = Blurred(loc(1):loc(2),loc(3):loc(4));
Blurred_gray_win = Blurred_gray(loc(1):loc(2),loc(3):loc(4));

figure;imshow(Blurred_gray);title('test Image');
figure;imshow(Blurred_gray_win);title('patch for kernel estimation');
%% iteration process
% build initial kernel for deblurring
blurstrength = 0.4;
initial_kernel = fspecial('gaussian', kernel_size , blurstrength);

delta = 0.0001;
err_image(1) = 0;

recent_dict = initdict;
recent_kernel = initial_kernel;

for i=1:iter_num
    if i>1
        previous_im = recent_image;
    end
    % update dictionary and recover image 
    param.blurry = Blurred_win;
    param.initdict = initdict;
    param.kernel = recent_kernel;
    param.rowskip = 3;
    param.colskip = 3;
    [recent_image, coef, recent_dict, lucy_image] = direct_sparse(param);
    figure;imshow(recent_image);title(['Iteration ',num2str(i)]);

    % estimate the kernel
    params.clear = recent_image;
    params.blurry = Blurred_win;
    params.lambda = 5;
    params.size = kernel_size;
    recent_kernel = find_kernel(params);
    
    if (i==1)
        err_image(i+1) = sqrt(sum(sum((recent_image-Blurred_gray_win).^2))/size(Blurred_win,1)/size(Blurred_win,2));
    else
        err_image(i+1) = sqrt(sum(sum((recent_image-previous_im).^2))/size(Blurred_win,1)/size(Blurred_win,2));
    end
    fprintf('Iteration: %i\n',i);
    fprintf('image error: %d\n',err_image(i+1));
    
    % break condition
    if (abs(err_image(i+1)-err_image(i)) < delta)
        break;
    end
end

% lucy algorithm to deblur the original image with estimated kernel
numit=20;
dampar=0.0005;
lim = ceil(size(recent_kernel,1)/2);
weight=zeros(size(ori_im));
weight(lim+1:end-lim,lim+1:end-lim)=1;

deblur_im = deconvlucy(Blurred,recent_kernel,numit,dampar,weight);
deblur_kernel = recent_kernel;