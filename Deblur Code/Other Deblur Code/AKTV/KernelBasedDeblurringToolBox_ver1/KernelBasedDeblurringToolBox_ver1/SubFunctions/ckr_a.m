function [z, zx, zy, zxx, zxy, zyy] = ckr_a(seq, mvs, r, ref, Q, h, ktype, ksize)
% Classic Kernel Regression
%
% [USAGE]
% [z zx zy] = ckr(seq, mvs, r, ref, Q, h, ktype, ksize)
% 
% [RETURNS]
% z      : an estimated image
% zx, zy : estimated image gradients along x and y directions
%
% [PARAMETERS]
% seq    : a video sequence
% mvs    : translational motion vectors
% ref    : a number of a reference frame
% r      : resolution enhancement factor
% Q      : regression order
%          0 = local constant estimator,
%          1 = local linear estimator,
%          2 = local quadratic estimator.
% h      : global smoothing parameter
% ktype  : type of the kernel function
%          'ga'  = gaussian,
%          'exp' = exponential (laplacian)
% ksize  : kernel size
%
% [HISTORY]
% Nov 19, 2005 : created by hiro
% Jan 20, 2006 : mirroring

% mirroring
[N M num] = size(seq);
width = floor(ksize / 2);
seqtmp = seq;
seq = 0;
seq = zeros(N + width * 2, M + width * 2, num);
for i = 1 : num
    seq(:,:,i) = EdgeMirror(seqtmp(:,:,i), [width, width]);
end
clear seqtmp;

a = 0;

% convert motion
[mvi, mvd] = convertmv(mvs, r, ref);

% dimension of low resolution frames
[N M num] = size(seq);
N = N * r;
M = M * r;

C = eye(2);

switch Q
    case 0 % Local Constant Estimator
        % initialize
        s11 = zeros(N, M);
        f1  = zeros(N, M);
        for n = 1 : num
            tmp = seq(:,:,n);
            if sum(tmp(:)) == 0
                continue;
            end
            I = upsample2(ones(size(seq(:,:,n))), r);
            shift = [mvd(2,n) mvd(1,n)];
            % obtain kernels
            K11 = lpkernel('11', ktype, ksize, h, C, a, shift);
            % densities
            Is = ishift2(I, mvi(2,n), mvi(1,n));
            s11 = s11 + conv2(Is, K11, 'same');
            % radiometric densities
            Yn = ishift2(upsample2(seq(:,:,n), r), mvi(2,n), mvi(1,n));
            f1 = f1 + conv2(Yn, K11, 'same');
        end
        zz = f1 ./ (s11 + (s11 == 0));
        z = zz(1+width*r : end-width*r, 1+width*r : end-width*r);
        zx = 0;
        zy = 0;
        return;
        
    case 1 % Local Linear Estimator
        global s11;
        global s12;
        global s22;
        global f1;
        global f2;
        % initialize
        s11 = zeros(N, M);
        s12 = zeros(N, M, 2);
        s22 = zeros(N, M, 4);
        f1  = zeros(N, M);
        f2  = zeros(N, M, 2);
        for n = 1 : num
            tmp = seq(:,:,n);
            if sum(tmp(:)) == 0
                continue;
            end
            I = upsample2(ones(size(seq(:,:,n))), r);
            shift = [mvd(2,n) mvd(1,n)];
            % obtain kernels
            K11 = lpkernel('11', ktype, ksize, h, C, a, shift);
            K12 = lpkernel('12', ktype, ksize, h, C, a, shift);
            K22 = lpkernel('22', ktype, ksize, h, C, a, shift);
            % densities
            Is = ishift2(I, mvi(2,n), mvi(1,n));
            s11 = s11 + conv2(Is, K11, 'same');
            s12(:,:,1) = s12(:,:,1) + conv2(Is, K12(:,:,1), 'same');
            s12(:,:,2) = s12(:,:,2) + conv2(Is, K12(:,:,2), 'same');
            s22(:,:,1) = s22(:,:,1) + conv2(Is, K22(:,:,1), 'same');
            s22(:,:,2) = s22(:,:,2) + conv2(Is, K22(:,:,2), 'same');
            s22(:,:,3) = s22(:,:,2);
            s22(:,:,4) = s22(:,:,4) + conv2(Is, K22(:,:,4), 'same');
            % radiometric densities
            Yn = ishift2(upsample2(seq(:,:,n), r), mvi(2,n), mvi(1,n));
            f1 = f1 + conv2(Yn, K11, 'same');
            f2(:,:,1) = f2(:,:,1) + conv2(Yn, K12(:,:,1), 'same');
            f2(:,:,2) = f2(:,:,2) + conv2(Yn, K12(:,:,2), 'same');
        end
        [zz zzx zzy] = beta_linear;
        zxx = 0;
        zxy = 0;
        zyy = 0;
        
    case 2 % Local Quadratic Estimator
        global s11;
        global s12;
        global s22;
        global s13;
        global s23;
        global s33;
        global f1;
        global f2;
        global f3;
        % initialize
        s11 = zeros(N, M);
        s12 = zeros(N, M, 2);
        s22 = zeros(N, M, 4);
        s13 = zeros(N, M, 3);
        s23 = zeros(N, M, 6);
        s33 = zeros(N, M, 9);
        f1  = zeros(N, M);
        f2  = zeros(N, M, 2);
        f3  = zeros(N, M, 3);
        for n = 1 : num
            tmp = seq(:,:,n);
            if sum(tmp(:)) == 0
                continue;
            end
            I = upsample2(ones(size(seq(:,:,n))), r);
            shift = [mvd(2,n) mvd(1,n)];
            % obtain kernels
            K11 = lpkernel('11', ktype, ksize, h, C, a, shift);
            K12 = lpkernel('12', ktype, ksize, h, C, a, shift);
            K22 = lpkernel('22', ktype, ksize, h, C, a, shift);
            K13 = lpkernel('13', ktype, ksize, h, C, a, shift);
            K23 = lpkernel('23', ktype, ksize, h, C, a, shift);
            K33 = lpkernel('33', ktype, ksize, h, C, a, shift);
            % densities
            Is = ishift2(I, mvi(2,n), mvi(1,n));
            s11 = s11 + conv2(Is, K11, 'same');
            s12(:,:,1) = s12(:,:,1) + conv2(Is, K12(:,:,1), 'same');
            s12(:,:,2) = s12(:,:,2) + conv2(Is, K12(:,:,2), 'same');
            s22(:,:,1) = s22(:,:,1) + conv2(Is, K22(:,:,1), 'same');
            s22(:,:,2) = s22(:,:,2) + conv2(Is, K22(:,:,2), 'same');
            s22(:,:,4) = s22(:,:,4) + conv2(Is, K22(:,:,4), 'same');
            s13(:,:,1) = s13(:,:,1) + conv2(Is, K13(:,:,1), 'same');
            s13(:,:,2) = s13(:,:,2) + conv2(Is, K13(:,:,2), 'same');
            s13(:,:,3) = s13(:,:,3) + conv2(Is, K13(:,:,3), 'same');
            s23(:,:,1) = s23(:,:,1) + conv2(Is, K23(:,:,1), 'same');
            s23(:,:,2) = s23(:,:,2) + conv2(Is, K23(:,:,2), 'same');
            s23(:,:,4) = s23(:,:,4) + conv2(Is, K23(:,:,4), 'same');
            s23(:,:,6) = s23(:,:,6) + conv2(Is, K23(:,:,6), 'same');
            s33(:,:,1) = s33(:,:,1) + conv2(Is, K33(:,:,1), 'same');
            s33(:,:,2) = s33(:,:,2) + conv2(Is, K33(:,:,2), 'same');
            s33(:,:,3) = s33(:,:,3) + conv2(Is, K33(:,:,3), 'same');
            s33(:,:,6) = s33(:,:,6) + conv2(Is, K33(:,:,6), 'same');
            s33(:,:,9) = s33(:,:,9) + conv2(Is, K33(:,:,9), 'same');
            % radiometric densities
            Yn = ishift2(upsample2(seq(:,:,n), r), mvi(2,n), mvi(1,n));
            f1 = f1 + conv2(Yn, K11, 'same');
            f2(:,:,1) = f2(:,:,1) + conv2(Yn, K12(:,:,1), 'same');
            f2(:,:,2) = f2(:,:,2) + conv2(Yn, K12(:,:,2), 'same');
            f3(:,:,1) = f3(:,:,1) + conv2(Yn, K13(:,:,1), 'same');
            f3(:,:,2) = f3(:,:,2) + conv2(Yn, K13(:,:,2), 'same');
            f3(:,:,3) = f3(:,:,3) + conv2(Yn, K13(:,:,3), 'same');
        end
        s22(:,:,3) = s22(:,:,2);
        s23(:,:,3) = s23(:,:,2);
        s23(:,:,5) = s23(:,:,4);
        s33(:,:,4) = s33(:,:,2);
        s33(:,:,5) = s33(:,:,3);
        s33(:,:,7) = s33(:,:,3);
        s33(:,:,8) = s33(:,:,6);
        
        [zz zzx zzy zzxx zzxy zzyy] = beta_quad2;
end

z = zz(1+width*r : end-width*r, 1+width*r : end-width*r);
zx = zzx(1+width*r : end-width*r, 1+width*r : end-width*r);
zy = zzy(1+width*r : end-width*r, 1+width*r : end-width*r);

if Q == 1
    return;
end

zxx = zzxx(1+width*r : end-width*r, 1+width*r : end-width*r);
zxy = zzxy(1+width*r : end-width*r, 1+width*r : end-width*r);
zyy = zzyy(1+width*r : end-width*r, 1+width*r : end-width*r);