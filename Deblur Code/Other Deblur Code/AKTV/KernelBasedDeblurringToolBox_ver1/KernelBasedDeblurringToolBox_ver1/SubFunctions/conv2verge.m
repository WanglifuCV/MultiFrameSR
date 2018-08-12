function z = conv2verge(img, kernel)
% CONV2VERGE [z = conv2verge(img, kernel)]
% this function computes 2 dimentional convolution of img and kernel with
% suppressing verge effect.
% kernel  : filter kernel
% img     : image
%
% coded by hiro on Nov 21, 2004

bksize = size(kernel);
imsize = size(img);
xm = round(bksize(1)/2);
ym = round(bksize(2)/2);

% add a reflection edge
tmp = EdgeMirror(img, [xm ym]);

tmp2 = conv2(tmp, kernel, 'same');
z = tmp2(xm+1:imsize(1)+xm, ym+1:imsize(2)+ym);