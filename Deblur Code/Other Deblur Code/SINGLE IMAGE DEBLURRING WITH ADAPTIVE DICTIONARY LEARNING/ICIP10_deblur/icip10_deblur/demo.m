addpath('ksvdbox') % K-SVD dictionary training algorithm
addpath('ompbox') % Orthogonal Matching Pursuit algorithm


ori_im = imread('blurred.jpg');
kernel_size = 5;
patch_size = 12;
iter_num = 8;
loc = [250,449,80,279];

[deblur_im,kernel] = deblur_adl(ori_im, kernel_size,patch_size,iter_num,loc);
figure;imshow(deblur_im);
