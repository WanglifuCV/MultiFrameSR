function K = lpkernel(order, ktype, size, h, C, a, shift)
% order : '11', '12', '13', '22', '23', or '33'

K = 0;
r = (size - 1) / 2.0;
[y, x] = meshgrid(-r:r, -r:r);
x = x - shift(1);
y = y - shift(2);
% x = x ./ 320;
% y = y ./ 320;
Cinv = inv(C * h^2);
%Cinv = (C / h^2);

tt = (x.*(Cinv(1,1)*x + Cinv(2,1)*y) + y.*(Cinv(1,2)*x + Cinv(2,2)*y)); % = t.^2

switch ktype
    case {'ep',  1} % Epanechnikov
        K = 1 - tt;
        K = K .* (K > 0);
    case {'bi',  2} % Biweight
        K = 1 - tt;
        K = K .* (K > 0);
        K = K.^2;
    case {'tri', 3} % Triangle
        K = 1 - sqrt(tt);
        K = K .* (K > 0);
    case {'exp', 4} % Exponential
        K = exp(-0.5 * sqrt(tt));
    case {'ga',  5} % Gaussian
        K = exp(-0.5 * tt);
    case {'bw',  6} % Butterworth
        K = 1 ./ (1 + tt.^a);
    case {'sk', 7} % spline kernel
        K = 0.5 .* exp(-sqrt(tt) / sqrt(2)) .* sin(sqrt(tt) / sqrt(2) + pi / 4);
    case {'tricube'} % tricube kernel
        K = 70 / 81 * (1 - sqrt(tt).^3).^3;
        K = K .* (sqrt(tt) <= 1);
end

% normalize
K = K / sqrt(det(C));

% %B = conv2(fspecial('ga', 11, 2.5), fspecial('ga', 5, 1.0), 'same');
% B = fspecial('ga', 5, 1.0);
% % x = conv2(x, B, 'same');
% % y = conv2(y, B, 'same');

tmp = K;
switch order
    case '11'
        K = tmp;
    case '12'
        K = zeros(size, size, 2);
        K(:,:,1) = x .* tmp;
        K(:,:,2) = y .* tmp;
    case '13'
        K = zeros(size, size, 3);
        K(:,:,1) = x.^2 .* tmp;
        K(:,:,2) = x .* y .* tmp;
        K(:,:,3) = y.^2 .* tmp;
    case '22'
        K = zeros(size, size, 4);
        K(:,:,1) = x.^2 .* tmp;
        K(:,:,2) = x .* y .* tmp;
        K(:,:,3) = K(:,:,2);
        K(:,:,4) = y.^2 .* tmp;
    case '23'
        K = zeros(size, size, 6);
        K(:,:,1) = x.^3 .* tmp;
        K(:,:,2) = x.^2 .* y .* tmp;
        K(:,:,3) = K(:,:,2);
        K(:,:,4) = x .* y.^2 .* tmp;
        K(:,:,5) = K(:,:,4);
        K(:,:,6) = y.^3 .* tmp;
    case '33'
        K = zeros(size, size, 9);
        K(:,:,1) = x.^4 .* tmp;
        K(:,:,2) = x.^3 .* y .* tmp;
        K(:,:,3) = x.^2 .* y.^2 .* tmp;
        K(:,:,4) = K(:,:,2);
        K(:,:,5) = K(:,:,3);
        K(:,:,6) = x .* y.^3 .* tmp;
        K(:,:,7) = K(:,:,3);
        K(:,:,8) = K(:,:,6);
        K(:,:,9) = y.^4 .* tmp;
end
% 
% tmp = K;
% switch order
%     case '11'
%         K = tmp;
%     case '12'
%         K = zeros(size, size, 2);
%         K(:,:,1) = conv2(x, B, 'same') .* tmp;
%         K(:,:,2) = conv2(y, B, 'same') .* tmp;
%     case '13'
%         K = zeros(size, size, 3);
%         K(:,:,1) = conv2(x.^2, B, 'same') .* tmp;
%         K(:,:,2) = conv2(x .* y, B, 'same') .* tmp;
%         K(:,:,3) = conv2(y.^2, B, 'same') .* tmp;
%     case '22'
%         K = zeros(size, size, 4);
%         K(:,:,1) = conv2(x, B, 'same') .* conv2(x, B, 'same') .* tmp;
%         K(:,:,2) = conv2(x, B, 'same') .* conv2(y, B, 'same') .* tmp;
%         K(:,:,3) = K(:,:,2);
%         K(:,:,4) = conv2(y, B, 'same') .* conv2(y, B, 'same') .* tmp;
%     case '23'
%         K = zeros(size, size, 6);
%         K(:,:,1) = conv2(x.^2, B, 'same') .* conv2(x, B, 'same') .* tmp;
%         K(:,:,2) = conv2(y, B, 'same') .* conv2(x .* x, B, 'same') .* tmp;
%         K(:,:,3) = conv2(x, B, 'same') .* conv2(x .* y, B, 'same') .* tmp;
%         K(:,:,4) = conv2(y, B, 'same') .* conv2(x .* y, B, 'same') .* tmp;
%         K(:,:,5) = conv2(x, B, 'same') .* conv2(y .* y, B, 'same') .* tmp;
%         K(:,:,6) = conv2(y, B, 'same') .* conv2(y.^2, B, 'same') .* tmp;
%     case '33'
%         K = zeros(size, size, 9);
%         K(:,:,1) = conv2(x.^2, B, 'same') .* conv2(x.^2, B, 'same') .* tmp;
%         K(:,:,2) = conv2(x .* y, B, 'same') .* conv2(x.^2, B, 'same') .* tmp;
%         K(:,:,3) = conv2(x.^2, B, 'same') .* conv2(y.^2, B, 'same') .* tmp;
%         K(:,:,4) = K(:,:,2);
%         K(:,:,5) = conv2(x .* y, B, 'same') .* conv2(x .* y, B, 'same') .* tmp;
%         K(:,:,6) = conv2(x .* y, B, 'same') .* conv2(y.^2, B, 'same') .* tmp;
%         K(:,:,7) = K(:,:,3);
%         K(:,:,8) = K(:,:,6);
%         K(:,:,9) = conv2(y.^2, B, 'same') .* conv2(y.^2, B, 'same') .* tmp;
% end


% make kernel size odd number
if mod(size, 2) == 0
    % if the kernel size is even numver, add a column and a row whose
    % elemens are all zero.
    [N M P] = size(K);
    tmp = zeros(N+1, M+1, P);
    tmp(1:N, 1:M, :) = K(:,:,:);
    K = 0;
    K = tmp;
end