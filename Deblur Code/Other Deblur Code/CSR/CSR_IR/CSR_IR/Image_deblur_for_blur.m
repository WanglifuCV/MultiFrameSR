% =========================================================================
% CSR-based image deblurring, Version 1.0
% Copyright(c) 2011 Weisheng Dong, Lei Zhang, Guangming Shi
% All Rights Reserved.
%
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is here
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%----------------------------------------------------------------------
%
% This is an implementation of the algorithm for non-blind image deblurring
% 
% Please refer to the following paper if you use this code:
%
% Weisheng Dong, Lei Zhang, and Guangming Shi,
% "Centralized sparse representation for image restoration", 
% in IEEE Int. Conf. Computer Vision (ICCV), 2011. 
% 
%--------------------------------------------------------------------------
clc;
clear;
addpath('Utilities');

blur_type          =    1;          % 1: uniform blur kernel;  2: Gaussian blur kernel; 
if blur_type == 1                   % When blur_type = 1, blur_par denotes the kernel size; When blur_type = 2, blur_par denotes the standard variance of Gaussian kernel
    blur_par       =    3;          % the default blur kernel size is 9 for uniform blur;
else
    blur_par       =    1.6;        % the default standard deviation of Gaussian blur kernel is 3
end
nSig               =    sqrt(2);    % The standard variance of the additive Gaussian noise;
out_dir            =    'Data\TestData\';
im_true_dir        =    'Data\TestData\';
im_true_name       =    'GroundTruth_0023.tif';
im_blur_dir        =    'Data\TestData\';
im_blur_name       =    'foreman_Proposed_SR_23.tif';


if blur_type==1
    if nSig<2
        par.t1        =   0.06*nSig^2;
        par.c1        =   0.07*nSig^2;    % 0.08*nSig^2;             
        par.t2        =   0.19*nSig^2;        
        par.c2        =   0.56*nSig^2;    %0.5*nSig^2;
        par.eps2      =   0.31;
    else
        par.t1        =   0.04*nSig^2;
        par.c1        =   0.06*nSig^2;
        par.t2        =   0.17*nSig^2;        
        par.c2        =   0.45*nSig^2;
        par.eps2      =   0.32;
    end
else
    par.t1        =   0.03*nSig^2;
    par.c1        =   0.04*nSig^2;        
    par.t2        =   0.12*nSig^2;        
    par.c2        =   0.36*nSig^2;
    par.eps2      =   0.31;
end
par.nSig      =   nSig;                                % Std variance of Gassian noise
par.iters     =   160;
par.nblk      =   12;
par.sigma     =   1.4;    
par.eps       =   0.3;

par.I                  =   double( imread(fullfile(im_true_dir, im_true_name)) );
bI                     =   double( imread(fullfile(im_blur_dir, im_blur_name)) );
figure(1); imshow(uint8(bI));
[par.bim par.fft_h]    =   Generate_blur_image_for_blur(bI, blur_type, blur_par, par.nSig);
[im PSNR SSIM]         =   CSR_Deblurring( par ); 

fname            =   strcat('CSR_', im_blur_name);    
imwrite(im./255, fullfile(out_dir, fname));    
disp( sprintf('%s: PSNR = %3.2f  SSIM = %f\n\n', im_blur_name, PSNR, SSIM) );


