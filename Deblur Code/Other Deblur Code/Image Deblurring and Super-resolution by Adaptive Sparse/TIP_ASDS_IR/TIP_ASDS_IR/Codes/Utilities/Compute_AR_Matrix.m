function  A     =  Compute_AR_Matrix( im, par )

Codeword   =  par.Codeword;
Modepar    =  par.AR_D;

cls_num    =  size(Codeword, 2);
f          =  3;       % The window size = f*2+1;
t          =  1;
ws         =  2*t+1;
cen        =  t+1;
tap_num    =  ws*ws-1;
t2         =  tap_num/2;

ext_im     =  padarray( im, [f f], 'symmetric' );
[h w]      =  size( ext_im );
LP_filter  =  fspecial('gaussian', 7,  par.sigma);
lp_im      =  conv2( LP_filter, ext_im );
hp_im      =  ext_im - lp_im(4:h+3, 4:w+3);


b     =  f*2+1;
N     =  h-b+1;
M     =  w-b+1;
L     =  N*M;
X     =  zeros(b*b,L);
k     =  0;
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
        blk  =  hp_im(i:h-b+i,j:w-b+j);
        blk  =  blk(:);
        X(k,:) =  blk';            
    end
end


[h, w]  =  size( im );
num_mod =  4;

pos     =  (1:h*w);
pos     =  reshape(pos, [h w]);

% Find the smooth blocks and set the smooth model
flag       =   ones(h*w, 1);
m          =   mean(X);
d          =   ( X - m( ones(size(X,1),1), :) ).^2;    % d          =   (X-repmat(m, size(X,1), 1)).^2;
v          =   sqrt( mean( d ) );
[a, ind]   =   find( v<par.delta );
flag(ind)  =   0;
flag       =   reshape(flag, h,w);
a          =   zeros(ws,ws);
a(cen-1:cen+1, cen-1:cen+1)  =  1;
a          =   a(:);
a(t2+1)    =   [];
a          =   a./sum(a);



% Set the triplet of the sparse matrix
nv      =  ws*ws;
nt      =  (nv)*h*w;
R       =  zeros(nt,1);
C       =  zeros(nt,1);
V       =  zeros(nt,1);
cnt     =  1;

for  row  = 1:h
    for  col  =  1:w
   
        off_cen =  (col-1)*h + row;
        
        if  flag(row,col)==0
            mode_par   =  a; 
        else        
            v          =  X(:,off_cen);
            dis        =  sum( (v(:,ones(1,cls_num) ) - Codeword).^2 );     %             wx         =  repmat( X(:, off_cen), 1, cls_num );            
            [m_dis Indics]  =  sort( dis );
                
            w_arr      =  1./( m_dis( Indics(1:num_mod) )+eps );
            w_arr      =  w_arr./sum( w_arr(:) );
        
            mode_par   =  repmat( w_arr, tap_num, 1 ).*Modepar(:, Indics(1:num_mod) );
            mode_par   =  sum( mode_par, 2 );                
        end
        
        rmin       =  max( row-t, 1);
        rmax       =  min( row+t, h);
        cmin       =  max( col-t, 1);
        cmax       =  min( col+t, w);
        sup        =  pos(rmin:rmax, cmin:cmax);
        col_ind    =  sup(:);
        
        r1         =  row-rmin;
        r2         =  rmax-row;
        c1         =  col-cmin;
        c2         =  cmax-col;
        
        AR         =  reshape([mode_par(1:t2); 0; mode_par(t2+1:end)], [ws, ws]);
        AR         =  -AR./sum(sum(AR(cen-r1:cen+r2, cen-c1:cen+c2)));
        AR(cen,cen)=  1;
        AR         =  AR(cen-r1:cen+r2, cen-c1:cen+c2);
        AR         =  AR(:);

        nn   = size(col_ind,1);
        
        R(cnt:cnt+nn-1)  =  off_cen;
        C(cnt:cnt+nn-1)  =  col_ind;
        V(cnt:cnt+nn-1)  =  AR;
        
        cnt              =  cnt + nn;
    end
end

R   =  R(1:cnt-1);
C   =  C(1:cnt-1);
V   =  V(1:cnt-1);
A   =  sparse(R, C, V, h*w, h*w);
