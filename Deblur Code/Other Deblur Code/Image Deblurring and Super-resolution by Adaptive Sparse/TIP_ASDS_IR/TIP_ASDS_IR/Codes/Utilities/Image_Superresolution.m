%----------------------------------------------------
% Image super-resolution
% Data: May 20th, 2010
% Author: Weisheng Dong, wsdong@mail.xidian.edu.cn
%----------------------------------------------------
function [im PSNR SSIM]   =  Image_Superresolution(method, dict, nSig, psf, scale, Output_dir, Test_image_dir, image_name)
par.method    =   method;                                % 0: ASDS;  1: ASDS_AR;  2: ASDS_AR_NL;
if dict==1
    load Data/Lib/PCA_AR_TD1;
    pre   =  'TD1_';
else    
    load Data/Lib/PCA_AR_TD2;
    pre   =  'TD2_';
end

if nSig == 0
    if par.method==0
        pre  =  strcat('ASDS_', pre);
        par.tau       =   0.08;    
        par.c1        =   0.7;
    elseif par.method==1
        pre  =  strcat('ASDS_AR_', pre);
        par.lam_AR    =   0.008;
        par.tau       =   0.08;
        par.c1        =   0.7;
    else
        pre  =  strcat('ASDS_AR_NL_', pre);
        par.lam_AR    =   0.008;
        par.lam_NL    =   0.04;
        par.tau       =   0.08;
        par.c1        =   0.7;
    end
    par.lamada    =   5.5;
else
    if par.method==0
        pre  =  strcat('ASDS_', pre);
        par.tau       =   0.7;
        par.c1        =   3.8;
    elseif par.method==1
        pre  =  strcat('ASDS_AR_', pre);
        par.lam_AR    =   0.06;
        par.tau       =   0.68;    
        par.c1        =   3.6;    
    else
        pre  =  strcat('ASDS_AR_NL_', pre);
        par.lam_AR    =   0.06;
        par.lam_NL    =   0.25;
        par.tau       =   0.66;
        par.c1        =   3.4;
    end
    par.lamada    =   0.8;
end

par.nSig      =   nSig;
par.scale     =   scale;
par.psf       =   psf;
par.nIter     =   999;

par.eps       =   2e-6;
par.PCA_D     =   PCA_D;
par.Codeword  =   centroids;
par.AR_D      =   AR_D;
par.nblk      =   10;

par.I         =    double( imread(fullfile(Test_image_dir, image_name)) );
LR            =    Blur('fwd', par.I, par.psf);
LR            =    LR(1:par.scale:end,1:par.scale:end,:);    
par.LR        =    Add_noise(LR, par.nSig);   
par.B         =    Set_blur_matrix( par );       
        
pp            =   'LR_';
fname         =    strcat(pp, image_name);
imwrite(par.LR./255, fullfile(Output_dir, fname));

pp           =   'NN_';
fname        =   strcat(pp, image_name);
NNim         =   imresize(par.LR, par.scale, 'nearest');
imwrite(NNim./255, fullfile(Output_dir, fname));
    

[im PSNR SSIM]   =   ASDS_AR_Superresolution( par );
fname            =   strcat(pre, image_name);
imwrite(im./255, fullfile(Output_dir, fname));
