function  im_out   =  LPCA_IS( im, PCA_idx, PCA_D, par )

b        =  par.win;
b2       =  b*b;
cls_num  =  size( PCA_D, 2 );
[h  w]   =  size(im);

N       =  h-b+1;
M       =  w-b+1;
L       =  N*M;
X       =  zeros(b*b, L);

% For the Y component
k    =  0;
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
        blk  =  im(i:end-b+i,j:end-b+j);
        blk  =  blk(:);
        X(k,:) =  blk';
    end
end

tau   =  par.tau;
Y     =  zeros( b2, L );
for i = 1:cls_num
       
    idx   =   find(PCA_idx==i);    
    if ( size(idx,1)==0 ) continue;  end
    
    P     =   reshape(PCA_D(:, i), b2, b2);
    xb    =   X(:, idx);
    yb    =   P'*soft(P*xb, tau);
    
    Y(:, idx)  =  yb;
end


% Output the processed image
im_out   =  zeros( size(im) );
im_wei   =  zeros(h,w);
k        =  0;
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
        blk  =  Y(k,:)';
        blk  =  reshape(blk, [N M]);
        im_out(i:h-b+i,j:w-b+j)  =  im_out(i:h-b+i,j:w-b+j) + blk;
        im_wei(i:h-b+i,j:w-b+j)  =  im_wei(i:h-b+i,j:w-b+j) + 1;
    end
end

im_out  =  im_out./(im_wei+eps);
