function  [PCA_D, D0, PCA_idx, s_idx, seg]  =  Set_PCA_idx( im, par, codewords, Arr )

[h  w]     =  size(im);
cls_num    =  size( codewords, 2 );
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
% X0      =  zeros(b*b, N*M, 'single');
X       =  zeros(b*b, L, 'single');
% Y       =  zeros(b*b, L, 'single');

% For the Y component
k    =  0;
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
%         blk  =  hp_im(i:end-b+i,j:end-b+j);  %blk  =  hp_im(r-1+i,c-1+j);
%         X0(k,:) =  blk(:)';
        
        blk    =  hp_im(r-1+i,c-1+j);
        X(k,:) =  blk(:)';        
    end
end
PCA_idx   =  zeros(L, 1);

m     =  mean(X);
d     =  ( X - m( ones(size(X,1),1), :) ).^2;     % d     =  (X-repmat(m, size(X,1), 1)).^2;
v     =  sqrt( mean( d ) );
[a, ind]  =   find( v<par.delta );
% D0        =   getpca( Y(:, ind) );

set         =  [1:L];
set(ind)    =  [];
L2          =  size(set,2);

for i = 1:L2
    
    wx            =   X(:, set(i));
    wx            =   wx(:, ones(1,cls_num));       %     wx            =   repmat( X(:,set(i)), 1, cls_num );
    dis           =   sum( (wx - codewords).^2 );        
    [md, idx]     =   min(dis);
    PCA_idx( set(i) )  =   idx;
%     dis    =  zeros(1, cls_num);
%     for j  =  1:8
%         wx            =   X0(:, Arr(j, set(i)));
%         wx            =   wx(:, ones(1,cls_num));       %     wx            =   repmat( X(:,set(i)), 1, cls_num );
%         dis               =   dis + sum( (wx - codewords).^2 );        
%     end
%     [md, idx]     =   min(dis);
%     PCA_idx( set(i) )  =   idx;

end

[s_idx, seg]  =   Proc_cls_idx( PCA_idx );
PCA_D         =   par.PCA_D(:, 2:end);
D0            =   reshape(par.PCA_D(:,1), b*b,b*b);     %  D0    =  dctmtx(b2);

% if flag==1
%     idx       =   s_idx(seg(1)+1:seg(2));
%     if  length(idx)>=200
%         D0        =   getpca( Y(:, idx) );
%     end
%     
%     for   i  = 2:length(seg)-1   
%         idx  =  s_idx(seg(i)+1:seg(i+1));    
%         if  length(idx) >= 180
%             cls  =  PCA_idx(idx(1));
%             pca  =   getpca( Y(:, idx) );
%             PCA_D(:, cls)   =  pca(:);
%         end
%     end
% end
return;
        