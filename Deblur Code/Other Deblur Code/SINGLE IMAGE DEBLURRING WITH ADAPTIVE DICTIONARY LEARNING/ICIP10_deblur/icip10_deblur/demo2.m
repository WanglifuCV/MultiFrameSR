function demo2()

ori_im = imread('castle.jpg');
kernel_size = 7;
patch_size = 12;
iter_num = 8;
loc = [40,209,70,239];

[deblur_im,kernel] = deblur_adl(ori_im, kernel_size,patch_size,iter_num,loc);
figure;imshow(deblur_im);
