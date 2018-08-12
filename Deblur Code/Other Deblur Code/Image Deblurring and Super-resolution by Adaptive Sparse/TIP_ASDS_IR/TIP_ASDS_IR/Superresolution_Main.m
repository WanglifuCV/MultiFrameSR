%----------------------------------------------------
% Sparse modeling with adaptive sparse domain selection 
% For Image Super-resolution
% Data: Nov. 20, 2010
% Author: Weisheng Dong, Lei Zhang, {wsdong@mail.xidian.edu.cn; cslzhang@comp.polyu.edu.hk}
%----------------------------------------------------
clc;
clear;
addpath('Codes');


Test_image_dir     =    'Data\SR_test_images';

psf                =     fspecial('gauss', 7, 1.6);              % The simulated PSF
scale              =    3;                                       % Downsampling factor 3
nSig               =    0;                                       % The standard variance of the additive Gaussian noise;
method             =    2;                                       % 0: ASDS;  1: ASDS_AR;  2: ASDS_AR_NL;
dict               =    1;                                       % 1: dictionary 1 trained from dataset 1; 2: dictionary 2 trained from dataset 2;

image_name         =    'Parrots.tif';                           % input the test image name;

if nSig == 0
    Output_dir         =    'Super_resolution\Noiseless_results';    % Where the output image will be generated;
else
    Output_dir         =    'Super_resolution\Noisy_results';        % Where the output image will be generated;
end

% The following codes start to perform the deblurring experiment with the above input parameters;

[im PSNR SSIM]   =   Image_Superresolution(method, dict, nSig, psf, scale, Output_dir, Test_image_dir, image_name);

disp( sprintf('%s: PSNR = %3.2f  SSIM = %f\n', image_name, PSNR, SSIM) ); 
