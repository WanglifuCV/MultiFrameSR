function fk = find_kernel(params, varargin)
% find_kernel function is to find the kernel using Tikhonov regularization
% method with recovered deblurred image and blurred image.
% Tikhonov regularization problem is a (convex) quadratic optimization problem: 
% minimize ||Ax-b||^2_2+\lambda ||x||_2^2 = x'(A'A+\lambda I)x-2b'Ax+b'b
% And the solution is:
% x= (A'A+\lambda I)^{-1}A'b
% 
% Required fields in PARAMS:
%  --------------------------
%    'clear' - the image we reconstructed from dictionary and coefficient
%    'blurry' - blurry image
%    'lambda' - the parameter to control the optimization problem
%    'size' - kernel size
% 
% Zhe Hu
% Computer Science Department
% University of California, Merced
% zhu@ucmerced.edu
% 
% October 2009

clear_image = params.clear;
blurry_image = params.blurry;
lambda = params.lambda;
size_k = params.size;

[clear_r, clear_c] =size(clear_image);
[blurry_r, blurry_c] =size(blurry_image);

if (clear_r ~= blurry_r) || (clear_c ~= blurry_c)
    error('recovered image and bllurred image have different size');
end

%% build A from recovered image
A = [];
for i=1:(clear_c-size_k+1)
    for j=1:(clear_r-size_k+1)
       patch = clear_image(j:j+size_k-1, i:i+size_k-1);
       patch = im2col(patch, [1 1], 'distinct');
       A = [A; patch];
    end
end

%% build b from blurry image
n = floor(size_k / 2);
b = blurry_image(n+1:blurry_r-n, n+1:blurry_c-n);
b = im2col(b, [1 1], 'distinct');
b = b';

%% find the solution
G = A' * A + lambda * eye(size_k * size_k);
if (det(G)==0)
G = inv (G' * G) * G';
else
G = inv(G);


fk = G * A' * b;

fk(isnan(fk)==1) = 10;
fk(fk<0) = 0;

% for i=1:size(fk,1)
%     if (fk(i)<0)
%         fk(i)=0;
%     end
% end

fk = fk / sum(fk(:));

fk = col2im(fk, [size_k size_k], [size_k size_k], 'distinct');
end