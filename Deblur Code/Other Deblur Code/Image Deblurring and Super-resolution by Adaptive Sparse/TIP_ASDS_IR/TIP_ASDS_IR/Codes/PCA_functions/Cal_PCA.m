function [PCA_D mX] =  Cal_PCA( im, par )

b    =  par.win;
t    =  floor( b/2 );
T    =  par.T;
st   =  par.step;
num  =  par.num;

[h  w]  =  size(im);

PCA_D   =  single(zeros(b*b*b*b, floor(h*w/st)));
mX      =  single(zeros(b*b, floor(h*w/st)));
cnt     =  1;

for  row  =  1:st:h

    if  (row+b-1) > h
        row   =  h-b+1;
    end       
    
    for  col  =  1:st:w
                
        if  (col+b-1) > w
            col   =  w-b+1;
        end                      

        cb      =  im(row:row+b-1, col:col+b-1);
           
        rmin    =  max( row+t-T,  1);
        rmax    =  min( row+t+T,  h);
        cmin    =  max( col+t-T,  1);
        cmax    =  min( col+t+T,  w);
        
        Blk     =  im(rmin:rmax, cmin:cmax);
            
        
        [n, m]  =  size(Blk);
        N       =  n-b+1;
        M       =  m-b+1;
        L       =  N*M;
        X       =  zeros(b*b,L);

        cb      =  cb';
        cb      =  cb(:);
        cb      =  repmat(cb,1,L);

        k   =  0;
        X   =  zeros(b*b,L);
        for i  = 1:b
            for j  = 1:b
                k    =  k+1;
                blk  =  Blk(i:n-b+i,j:m-b+j);
                blk  =  blk(:);
                X(k,:) =  blk';
            end
        end


        E    =  abs(X-cb).^2;
        mE   =  mean(E);
        [val, ind]   =  sort(mE);
        
        num  =  min( par.num, L );
        X    =  X(:,ind(1:num));
        
        [P, mx]   =  getpca(X);
        
        PCA_D(:,cnt)   =  single(P(:));
        mX(:,cnt)      =  single(mx);
        cnt            =  cnt + 1;
                
    end
    
end
        