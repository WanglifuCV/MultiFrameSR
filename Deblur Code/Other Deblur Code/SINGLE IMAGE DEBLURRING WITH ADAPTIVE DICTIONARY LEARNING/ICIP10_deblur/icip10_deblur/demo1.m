addpath('ksvdbox13') % K-SVD dictionary training algorithm
addpath('ompbox10') % Orthogonal Matching Pursuit algorithm

ori_im = imread('koala.jpg');
kernel_size = 7;
patch_size = 12;
iter_num = 8;
loc = [30,179,80,229];

[deblur_im,kernel] = deblur_adl(ori_im, kernel_size,patch_size,iter_num,loc);
figure;imshow(deblur_im);