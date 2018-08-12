% =========================================================================
% CSR-based image restoration, Version 1.0
% Copyright(c) 2011 Weisheng Dong, Lei Zhang, Guangming Shi
% All Rights Reserved.
%
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is here
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%----------------------------------------------------------------------
%
% This is an implementation of the CSR regularization algorithm
% 
% Please refer to the following paper if you use this code:
%
% Weisheng Dong, Lei Zhang, Guangming Shi,"Centralized sparse
% representation for image restoration", to appear in IEEE Int. Conf.
% Computer Vision (ICCV), 2011. 
% 
%--------------------------------------------------------------------------
function  [im_out]  =  CSR_Regularization( im, opts, Dict, blk_arr, wei_arr, lam, gam, flag )
[h w ch]   =   size(im);
b          =   opts.win;
b2         =   b*b*ch;
k          =   0;
s          =   opts.step;
A          =   Dict.D0;
PCA_D      =   Dict.PCA_D;
PCA_idx    =   Dict.cls_idx;
s_idx      =   Dict.s_idx;
seg        =   Dict.seg;

N     =  h-b+1;
M     =  w-b+1;
L     =  N*M;
r     =  [1:s:N];
r     =  [r r(end)+1:N];
c     =  [1:s:M];
c     =  [c c(end)+1:M];
X     =  zeros(b*b,L,'single');
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
        blk  =  im(i:h-b+i,j:w-b+j);
        blk  =  blk(:);
        X(k,:) =  blk';            
    end
end

X_m     =  zeros(length(r)*length(c),b*b,'single');
X0 = X';
for i = 1:opts.nblk
   v        =  wei_arr(:,i);
   X_m      =  X_m(:,:) + X0(blk_arr(:,i),:) .*v(:, ones(1,b2));
end
X_m=X_m';

ind        =   zeros(N,M);
ind(r,c)   =   1;
X          =   X(:, ind~=0);
N          =   length(r);
M          =   length(c);
L          =   N*M;
Y          =   zeros( b2, L );

idx        =   s_idx(seg(1)+1:seg(2));
t1         =   opts.t1;
t2         =   opts.t2;
if  flag==1
    t1    =   lam(:, idx);
    t2    =   gam(:, idx);
end
Y(:, idx)    =   A'*Biv_thresholding(A*X(:, idx), A*X_m(:, idx), t1, t2); 

for   i  = 2:length(seg)-1   
    idx    =   s_idx(seg(i)+1:seg(i+1));    
    cls    =   PCA_idx(idx(1));
    P      =   reshape(PCA_D(:, cls), b2, b2);
    t1     =   opts.t1;
    t2     =   opts.t2;
    if  flag==1 
        t1    =   lam(:, idx);
        t2    =   gam(:, idx);        
    end
    Y(:, idx)  =  P'*Biv_thresholding(P*X(:, idx), P*X_m(:, idx), t1, t2);
end

im_out   =  zeros(h,w);
im_wei   =  zeros(h,w);
k        =  0;
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
        im_out(r-1+i,c-1+j)  =  im_out(r-1+i,c-1+j) + reshape( Y(k,:)', [N M]);
        im_wei(r-1+i,c-1+j)  =  im_wei(r-1+i,c-1+j) + 1;       
    end
end
im_out    =  im_out./(im_wei+eps);
return;


function Y  = Biv_thresholding( x, mu, t1, t2 )
m1  =  (mu>=0).*mu;
m2  =  (mu<0).*mu;

y1   =  ( x < (-t1-t2) ).*( x+t1+t2 );
y3   =  ( x>(t1-t2) & x<(t1-t2+m1) ).*( x-t1+t2 );
y4   =  ( x>=(t1-t2+m1) & x<=(t1+t2+m1) ).*m1;
y5   =  ( x>(t1+t2+m1) ).*( x-t1-t2 );
Y1   =  (y1 + y3 + y4 + y5).*(mu>=0);

y1   =  ( x<(m2-t1-t2) ).*( x+t1+t2 );
y2   =  ( x>=(m2-t1-t2) & x<=(t2+m2-t1) ).*m2;
y3   =  ( x>(t2+m2-t1) & x<(t2-t1) ).*( x+t1-t2 );
y5   =  ( x>(t1+t2) ).*( x-t1-t2 );
Y2   =  (y1 + y2 + y3 + y5).*(mu<0);
Y    =  Y1 + Y2;
return;
