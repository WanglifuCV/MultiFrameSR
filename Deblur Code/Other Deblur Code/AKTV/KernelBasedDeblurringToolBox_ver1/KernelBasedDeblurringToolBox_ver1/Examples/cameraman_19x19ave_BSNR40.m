%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEERING KERNEL DEBLURRING EXAMPLE 1 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This simulation example demonstrates the kernel regression-based
% deblurring method with the steering kernel function, and this generates
% the deblurred image of Fig.10(g) in the paper "Deblurring Using
% Locally-Adaptive Kernel regression.
%
% [Details]
% test image : cameraman
% PSF : 19 x 19 uniform
% BSNR : 40[dB]
%
% [History]
% Oct 12, 2007 : coded and debugged by Hiro

% the mirroring width
mg = 15;

% read a test image
img = double(imread('cameraman.tif')); % Cameraman image

% create a PSF
PSFsupport = 19;
A = fspecial('ave', PSFsupport); % 19 x 19 uniform

% blur the test image
imgb = conv2verge(img, A);

% add white Gaussian noise
BSNR = 40; % blurred signal to noise ration in [dB]
vn = var(imgb(:)) / 10^(BSNR/10); % compute the noise variance
randn('state', 0); % initialize the noise generator 
y = imgb + randn(size(imgb)) * sqrt(vn); % add white noise to the blurred image

% mirror the edges
y = EdgeMirror(y, [mg, mg]);
[N M] = size(y);

% initialization
% the pilot estimation by weiner filter
K = 0.0001;
Af = fft2(A, 512, 512);
F = (abs(Af).^2 ./ (abs(Af).^2 + K)) ./ Af .* fft2(EdgeMirror(y(1+mg:end-mg, 1+mg:end-mg), [128, 128]));
y_init = ishift2(real(ifft2(F)), 9, 9);
y_init = y_init(129-mg:end-128+mg, 129-mg:end-128+mg);
y_init = double(round0_255(y_init));

fig_h = figure;
imshow(uint8(y_init(1+mg:end-mg, 1+mg:end-mg))); colormap(gray); axis image;
pause(0.1);

% initialization by the second order classic kernel regression
title('Initializing by the second order classic kernel regression...');
pause(0.1);
Q = 2; % the regression order
hc = 0.8; % the global smoothing parameter
ktype = 'ga'; % the kernel type ('ga' = Gaussian)
ksize = 7; % the kernel support size
[U Ux Uy Uxx Uxy Uyy] = ckr_a(y_init, [0;0], 1, 1, Q, hc, ktype, ksize);
U = y_init;

% parameter
h = 0.5; % the global smoothing parameter for the likelihood term
hu = 1.2; % the global smoothing parameter for the regularization term
nu = 0.4; % the step size
IT = 6500; % the number of iterations
ksize = 7; % the kernel support size
ktype = 'ga'; % the kernel support size
reg_lambda = 0.75; % the regularization parameter
radius = (ksize - 1) / 2;

% parameters for the steering matrices
wsize = 11; % the analysis window size
lambda = 1; % the regularization parameter
alpha = 0.5; % the structure sensitivity
r = 1;

[x2, x1] = meshgrid(-radius:radius, -radius:radius);
x1sq = x1.^2;
x2sq = x2.^2;
x1x2 = x1.* x2;

for it = 1 : IT
    % initialize error corrections
    ec_U = zeros(N, M);
    ec_Ux = zeros(N, M);
    ec_Uy = zeros(N, M);
    
    % every 50 iterations, we re-create the weight matrices
    if mod(it-1, 50) == 0
        % estimate the orientation information from the estimated gradients
        % for the regularization term
        title('Estimating the orientation information from the estimated gradients...');
        pause(0.1);
        C = steering_modified(Ux, Uy, wsize, lambda, alpha, 1);
        % estimate the orientation information from the estimated blurred
        % gradients for the likelihood term
        title('Estimating the orientation information from the estimated blurred gradients...');
        pause(0.1);
        Cb = steering_modified(conv2verge(Ux, A), conv2verge(Uy, A), wsize, lambda, alpha, 1);
        C11 = zeros(N, M);
        C12 = zeros(N, M);
        C22 = zeros(N, M);
        detC = zeros(N, M);
        Cb11 = zeros(N, M);
        Cb12 = zeros(N, M);
        Cb22 = zeros(N, M);
        detCb = zeros(N, M);
        title('Creating the weight matrices...');
        pause(0.1);
        for n = 1 : N
            for m = 1 : M
                Ctmp = inv(C(:,:,n,m)) ./ hu^2;
                C11(n,m) = Ctmp(1,1);
                C12(n,m) = Ctmp(1,2);
                C22(n,m) = Ctmp(2,2);
                detC(n,m) = sqrt(det(Ctmp));
                Ctmp = inv(Cb(:,:,n,m)) ./ h^2;
                Cb11(n,m) = Ctmp(1,1);
                Cb12(n,m) = Ctmp(1,2);
                Cb22(n,m) = Ctmp(2,2);
                detCb(n,m) = sqrt(det(Ctmp));
            end
        end
        clear C;
        clear Cb;
        C11 = EdgeMirror(C11, [radius, radius]);
        C12 = EdgeMirror(C12, [radius, radius]);
        C22 = EdgeMirror(C22, [radius, radius]);
        detC = EdgeMirror(detC, [radius, radius]);
        Cb11 = EdgeMirror(Cb11, [radius, radius]);
        Cb12 = EdgeMirror(Cb12, [radius, radius]);
        Cb22 = EdgeMirror(Cb22, [radius, radius]);
        detCb = EdgeMirror(detCb, [radius, radius]);
        % compute weights
        Wu = zeros(N * r, M * r, ksize, ksize);
        Wz = zeros(N * r, M * r, ksize, ksize);
        for n = 1 : N * r
            for m = 1 : M * r
                  xCx =  x1sq .* C11(n:n+ksize-1, m:m+ksize-1) ...
                       + x2sq .* C22(n:n+ksize-1, m:m+ksize-1) ...
                       + x1x2 .* C12(n:n+ksize-1, m:m+ksize-1) .* 2;
                       
                  W = exp(-0.5 * xCx) .* detC(n:n+ksize-1, m:m+ksize-1);
                  Wu(n,m,:,:) = reshape(W, 1, 1, ksize, ksize);
                  xCx =  x1sq .* Cb11(n:n+ksize-1, m:m+ksize-1) ...
                       + x2sq .* Cb22(n:n+ksize-1, m:m+ksize-1) ...
                       + x1x2 .* Cb12(n:n+ksize-1, m:m+ksize-1) .* 2;    
                  W = exp(-0.5 * xCx) .* detCb(n:n+ksize-1, m:m+ksize-1);
                  Wz(n,m,:,:) = reshape(W, 1, 1, ksize, ksize);
            end
        end
        clear C11;
        clear C12;
        clear C22;
        clear detC;
        clear Cb11;
        clear Cb12;
        clear Cb22;
        clear detCb;
        Wu = Wu ./ max([max(Wu(:)); max(Wz(:))]);
        Wz = Wz ./ max([max(Wu(:)); max(Wz(:))]);
    end
    
    % compute the gradient term (dC(U)/dU) of the steepest descent
    for i = -radius : 1 : radius
        for j = -radius : 1 : radius
            
            dif = zeros(N * r, M * r);

            % gradients of the likelihood term
            Utmp = U + Ux*i + Uy*j;
            wt = Wz(:,:,i+radius+1,j+radius+1); % obtain the weights
            tmp = (y - conv2(ishift2(Utmp, -i, -j), A, 'same')) .* wt; % (Y - SSGU)
            Y_SSGU = zeros(N * r, M * r);
            Y_SSGU(1+mg:end-mg, 1+mg:end-mg) = tmp(1+mg:end-mg, 1+mg:end-mg);
            dif = dif + conv2(ishift2(Y_SSGU, i, j), A, 'same');
            dif2 = dif;
            
            % gradients of the regularization term
            wt = Wu(:,:,i+radius+1,j+radius+1) .* reg_lambda;
            U_SSIU = sign(U - ishift2(Utmp, -i, -j));
            WU_SSIU = U_SSIU .* wt;
            SSWU_SSIU = ishift2(WU_SSIU, i, j);
            dif = dif - (WU_SSIU - SSWU_SSIU);
            dif2 = dif2 + SSWU_SSIU;
            
            tmp = dif;
            dif = zeros(N * r, M * r);
            dif(1+mg:end-mg, 1+mg:end-mg) = tmp(1+mg:end-mg, 1+mg:end-mg);
            tmp = dif2;
            dif2 = zeros(N * r, M * r);
            dif2(1+mg:end-mg, 1+mg:end-mg) = tmp(1+mg:end-mg, 1+mg:end-mg);
            
            ec_U = ec_U + dif;
            ec_Ux = ec_Ux + dif2 * i;
            ec_Uy = ec_Uy + dif2 * j;
        end
    end

    ec_U_max = max(abs(ec_U(:)));
    step = nu / ec_U_max;    
    
    % update the image and its gradients 
    U = U + step * ec_U;
    Ux = Ux + step * ec_Ux;
    Uy = Uy + step * ec_Uy;
    U = EdgeMirror(U(1+mg:end-mg, 1+mg:end-mg), [mg, mg]);
    Ux = EdgeMirror(Ux(1+mg:end-mg, 1+mg:end-mg), [mg, mg]);
    Uy = EdgeMirror(Uy(1+mg:end-mg, 1+mg:end-mg), [mg, mg]);
    
    % compute the root mean square error
    error = img - U(1+mg:N-mg, 1+mg:M-mg);
    rmse = sqrt(mse(error));
    
    % display the image
    figure(fig_h); imshow(uint8(U(1+mg:N-mg, 1+mg:M-mg))); colormap(gray); axis image;
    title(['IT = ', num2str(it), ', RMSE = ', num2str(rmse)]);
    pause(0.1);
end

save cameraman_19x19ave_BSNR40_L2L1_N1N1_RMSE140278 U Ux Uy;

