clear; clc; close all;
% y = im2double(imread('result_SR_1.tif'));
y = im2double(imread('temp.tif'));
yn   =  y + 5/255*randn(size(y));
figure; imshow(yn);
[NA, y_est] = BM3D(1, yn, 10);
figure; imshow(y);   
figure; imshow(y_est);

% % y = im2double(imread('result_SR_1.tif'));
% y = im2double(imread('temp.tif'));
% % [NA, y_est] = BM3D(1, y, 10);
% figure; imshow(y);   
% % figure; imshow(y_est);
% I = im2double(imread('foreman_Proposed_SR_23_dnoise.tif'));
% figure; imshow(I);
% resd = (I-y).^2;
% P = resd > max(resd(:)) * 10^-2;
% II = y;
% II(P) = I(P);
% figure; imshow(II);
% 
% 
% % imwrite(y_est, 'temp_dnoise.tif');

