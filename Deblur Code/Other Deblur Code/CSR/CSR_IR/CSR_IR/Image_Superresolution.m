% =========================================================================
% CSR-based image super-resolution, Version 1.0
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
% This is an implementation of the algorithm for image super-resolution
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

psf                =    fspecial('gauss', 7, 1.6);
nSig               =    0;
scale_factor       =    3;
out_dir            =    'Results\SR_results\';
im_dir             =    'Data\SR_test_images\';
im_name            =    'Butterfly.tif';

if nSig==0
    par.t1        =   0.08;    
    par.t2        =   0.3; 
    par.c1        =   0.2;
    par.c2        =   1.4;
    par.lamada    =   7;
    par.n         =   5;
else
    par.t1        =   0.002*nSig^2;    
    par.t2        =   0.017*nSig^2; 
    par.c1        =   0.03;
    par.c2        =   0.19;
    par.lamada    =   1.0;
    par.n         =   3;
end
par.psf       =   psf;
par.scale     =   scale_factor;
par.nSig      =   nSig;
par.iters     =   160;
par.nblk      =   12;
par.sigma     =   1.4;    
par.eps       =   0.3;

par.I        =   double( imread(fullfile(im_dir, im_name)) );
LR           =   Blur('fwd', par.I, par.psf);
LR           =   LR(1:par.scale:end,1:par.scale:end,:);    
par.LR       =   Add_noise(LR, par.nSig);   
par.B        =   Set_blur_matrix( par );

[im PSNR SSIM]   =   CSR_Superresolution( par );  

fname            =   strcat('CSR_', im_name);    
imwrite(im./255, fullfile(out_dir, fname));    
disp( sprintf('%s: PSNR = %3.2f  SSIM = %f\n\n', im_name, PSNR, SSIM) );


