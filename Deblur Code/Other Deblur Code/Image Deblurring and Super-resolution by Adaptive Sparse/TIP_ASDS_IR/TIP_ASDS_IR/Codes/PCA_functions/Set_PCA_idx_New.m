function  [PCA_idx, s_idx, seg]  =  Set_PCA_idx_New( im, par, PCA_D, DCT )

[h  w]     =  size(im);
cls_num    =  size( PCA_D, 2 );
LP_filter  =  fspecial('gaussian', 7, par.sigma);
lp_im      =  conv2( LP_filter, im );
lp_im      =  lp_im(4:h+3, 4:w+3);
hp_im      =  im - lp_im;


b       =  par.win;
s       =  par.step;
N       =  h-b+1;
M       =  w-b+1;

r       =  [1:s:N];
r       =  [r r(end)+1:N];
c       =  [1:s:M];
c       =  [c c(end)+1:M];
L       =  length(r)*length(c);
X       =  zeros(b*b, L, 'single');
Y       =  zeros(b*b, L, 'single');

% For the Y component
k    =  0;
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
        blk  =  hp_im(r-1+i,c-1+j);
        Y(k,:) =  blk(:)';
        
        blk  =  im(r-1+i,c-1+j);
        X(k,:) =  blk(:)';        
    end
end
PCA_idx   =  zeros(L, 1);

m     =  mean(Y);
d     =  (Y-repmat(m, size(Y,1), 1)).^2;
v     =  sqrt( mean( d ) );
[a, ind]  =  find( v<par.delta );
set         =  [1:L];
set(ind)    =  [];
L2          =  size(set,2);


m           =  ceil(size(X,1)/2);
b2          =  b*b;
E           =  zeros(cls_num, L2, 'single');

% Y           =  DCT(1:m,:)*X;
% E(1,:)      =  sum(Y.^2);

for i = 1:cls_num
    
    P        =  reshape(PCA_D(:, i), b2, b2);
    Y        =  P(1:m,:)*X(:,set);
    E(i,:)   =  sum(Y.^2);    
end

[a, ind]   =  max(E);
PCA_idx(set) = ind';
% PCA_idx    =  (ind-1)';

[s_idx, seg] =   Proc_cls_idx( PCA_idx );
        