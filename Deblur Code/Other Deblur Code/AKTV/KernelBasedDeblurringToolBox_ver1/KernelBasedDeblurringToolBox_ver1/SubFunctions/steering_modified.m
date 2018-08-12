function C = steering_modified(zx, zy, wsize, lambda, alpha, r)
% Compute steering matrices
%
% [USAGE]
% C = steering(zx, zy, wsize, lambda, alpha)
%
% [RETURNS]
% C      : steering matrices
%
% [PARAMETERS]
% zx, zy : image gradients along x and y directions
% wsize  : size of an analysis window
% lambda : regularization parameter
% alpha  : structure sensitive parameter
%
% [HISTORY]
% Nov 19, 2005 : created by hiro
% June 15, 2007 : added spatial kernel

[N M] = size(zx);
C = zeros(2, 2, N, M);

if mod(wsize, 2) == 0
    wsize = wsize + 1;
end
win = (wsize - 1) / 2;

% spatial weight kernel
K = fspecial('disk', win);
K = K ./ max(K(:));

% mirroring
zx = EdgeMirror(zx, [win, win]);
zy = EdgeMirror(zy, [win, win]);

for i = 1 : r : N
    ii = (i - 1) / r + 1;
    for j = 1 : r : M
        jj = (j - 1) / r + 1;
        gx = 0;
        gy = 0;
        gx = zx(i:i+win*2, j:j+win*2) .* K;
        gy = zy(i:i+win*2, j:j+win*2) .* K;

        G = [gx(:), gy(:)];
        %len = length(gx(:));
        len = sum(K(:));
        [u s v] = svd(G,0);
        S(1) = ((s(1,1) + lambda) / (s(2,2) + lambda));
        S(2) = 1/S(1);
        tmp = (S(1) * v(:,1) * v(:,1).' + S(2) * v(:,2) * v(:,2).')  * ((s(1,1) * s(2,2) + 0.0000001) / len)^alpha;

        C(1:2, 1:2, ii, jj) = inv(tmp + (eye(2)*0.00001));
    end
end