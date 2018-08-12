%----------------------------------------------------
% ASDS_X image deblurring
% Data: May 20th, 2010
% Author: Weisheng Dong, wsdong@mail.xidian.edu.cn
%----------------------------------------------------
function  [im PSNR SSIM]   =  Image_Deblurring(method, dict, nsig, Out_dir, Test_image_dir, image_name, blur_type, blur_par, blur_dir, blur_name)

par.method    =   method;                                % 0: ASDS;  1: ASDS_AR;  2: ASDS_AR_NL;
if dict==1
    load Data/Lib/PCA_AR_TD1;
    pre   =  'TD1_';
elseif dict==2    
    load Data/Lib/PCA_AR_TD2;
    pre   =  'TD2_';
end

if blur_type==1                             % using different parameter setting for different blur kernel
    if par.method==0
        pre  =  strcat('ASDS_', pre);
        par.tau       =   0.07*nsig^2;    %0.12*nsig^2;
        par.c1        =   0.64*nsig^2;
    elseif par.method==1
        pre  =  strcat('ASDS_AR_', pre);
        par.gam       =   0.003;
        par.tau       =   0.07*nsig^2;    %0.1*nsig^2;
        par.c1        =   0.64*nsig^2;
    else
        pre  =  strcat('ASDS_AR_NL_', pre);
        par.gam       =   0.002;
        par.eta       =   0.018;
        par.tau       =   0.07*nsig^2;    
        par.c1        =   0.64*nsig^2;
    end
else
    if par.method==0
        pre  =  strcat('ASDS_', pre);
        par.tau       =   0.07*nsig^2;
        par.c1        =   0.4*nsig^2;
    elseif par.method==1
        pre  =  strcat('ASDS_AR_', pre);
        par.gam       =   0.006;
        par.tau       =   0.05*nsig^2; 
        par.c1        =   0.4*nsig^2;
    else
        pre  =  strcat('ASDS_AR_NL_', pre);
        par.gam       =   0.002;
        par.eta       =   0.018;
        par.tau       =   0.028*nsig^2;    
        par.c1        =   0.22*nsig^2;    
    end    
end

par.PCA_D     =   PCA_D;
par.Codeword  =   centroids;
par.AR_D      =   AR_D;

par.nSig      =   nsig;
par.nIter     =   999;
par.eps       =   1e-8;
par.nblk      =   15;
par.delta     =   4;

par.I         =    double( imread(fullfile(Test_image_dir, image_name)) );
bI            =    double( imread(fullfile(blur_dir, blur_name)) );
figure(1); imshow(uint8(bI));
[par.bim par.fft_h]     =   Generate_blur_image(bI, blur_type, blur_par, par.nSig);
    
[im PSNR SSIM]   =   ASDS_AR_deblurring( par );        
fname            =   strcat(pre, image_name);
imwrite(im./255, fullfile(Out_dir, fname));
figure(2); imshow(uint8(im));

