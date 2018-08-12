function  [Py  Px]  =  Get_patches( im, b, psf )

[h w ch]  =  size(im);
ws        =  floor( size(psf,1)/2 );

if  ch==3
    lrim      =  rgb2ycbcr( uint8(im) );
    im        =  double( lrim(:,:,1));    
end

lp_im     =  conv2( psf, im );
lp_im     =  lp_im(ws+1:h+ws, ws+1:w+ws);
hp_im     =  im - lp_im;

N         =  h-b+1;
M         =  w-b+1;
s         =  1;
r         =  [1:s:N];
r         =  [r r(end)+1:N];
c         =  [1:s:M];
c         =  [c c(end)+1:M];
L         =  length(r)*length(c);
Py        =  zeros(b*b, L, 'single');
Px        =  zeros(b*b, L, 'single');

k    =  0;
for i  = 1:b
    for j  = 1:b
        k       =  k+1;
        blk     =  hp_im(r-1+i,c-1+j);
        Py(k,:) =  blk(:)';
        
        blk     =  im(r-1+i,c-1+j);
        Px(k,:) =  blk(:)';        
    end
end

rand('seed',0);
P    =   randperm(L)';
N    =   round(L*0.6);
P    =   P(1:N);
Py   =   Py(:, P);
Px   =   Px(:, P);