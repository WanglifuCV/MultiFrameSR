function  im_out   =  LPCA_IS_fast( im, PCA_idx, PCA_D, par, s_idx, seg, A, flag, Tau, Arr, Wei )
[h w ch]   =   size(im);

b          =   par.win;
b2         =   b*b;
k          =   0;
s          =   par.step;

N     =  h-b+1;
M     =  w-b+1;
r     =  [1:s:N];
r     =  [r r(end)+1:N];
c     =  [1:s:M];
c     =  [c c(end)+1:M];
X     =  zeros(b*b,N*M,'single');
X_m   =  zeros(b*b,length(r)*length(c),'single');
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
        blk  =  im(i:h-b+i,j:w-b+j);
        X(k,:) =  blk(:)';            
    end
end

% Compute the mean blks
idx      =   s_idx(seg(1)+1:seg(2));
set      =   1:size(X_m,2);
set(idx) =   [];

for i = 1:par.nblk
   v            =  Wei(i,set);
   X_m(:,set)   =  X_m(:,set) + X(:, Arr(i, set)) .*v(ones(b2,1), :);
end


% X_m     =  X_m';
% X0 = X';
% for i = 1:par.nblk
%    v            =  Wei(set,i);
%    X_m(set,:)   =  X_m(set,:) + X0(Arr(set,i),:) .*v(:, ones(1,b2));
% end
% X_m=X_m';


ind         =   zeros(N,M);
ind(r,c)    =   1;
X           =   X(:, ind~=0);

N           =   length(r);
M           =   length(c);
L           =   N*M;
Y           =   zeros( b2, L );


% Smooth blocks

tau         =   par.tau;
if  flag==1
   tau    =  Tau(:, idx); 
end
Y(:, idx)    =   A'*soft( A*X(:,idx), tau );

for   i  = 2:length(seg)-1
   
    idx  =  s_idx(seg(i)+1:seg(i+1));
    
    cls  =  PCA_idx(idx(1));
    P    =   reshape(PCA_D(:, cls), b2, b2);

    tau   =  par.tau;
    if  flag==1
       tau    =  Tau(:, idx); 
    end
    Y(:, idx)    =   P'*soft(P*(X(:, idx)-X_m(:,idx)), tau) + X_m(:,idx);
end


% Output the processed image
im_out   =  zeros(h,w);
im_wei   =  zeros(h,w);
k        =  0;
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
        im_out(r-1+i,c-1+j)  =  im_out(r-1+i,c-1+j) + reshape( Y(k,:)', [N M]);
        im_wei(r-1+i,c-1+j)  =  im_wei(r-1+i,c-1+j) + 1;       
        
%         im_out(i:h-b+i,j:w-b+j)  =  im_out(i:h-b+i,j:w-b+j) + reshape( Y(k,:)', [N M]);
%         im_wei(i:h-b+i,j:w-b+j)  =  im_wei(i:h-b+i,j:w-b+j) + 1;
    end
end

im_out  =  im_out./(im_wei+eps);



% function  im_out   =  LPCA_IS_fast( im, PCA_idx, PCA_D, par, s_idx, seg, A, flag, Tau )
% 
% b        =  par.win;
% s        =  par.step;
% b2       =  b*b;
% [h  w]   =  size(im);
% 
% N       =  h-b+1;
% M       =  w-b+1;
% r       =  [1:s:N];
% r       =  [r r(end)+1:N];
% c       =  [1:s:M];
% c       =  [c c(end)+1:M];
% 
% N       =  length(r);
% M       =  length(c);
% X       =  zeros(b2, N*M);
% Y       =  X;
% 
% % For the Y component
% k    =  0;
% for i  = 1:b
%     for j  = 1:b
%         k    =  k+1;
%         blk  =  im(r-1+i,c-1+j);
%         X(k,:) =  blk(:)';
%     end
% end
% 
% e1          =   0.2;
% 
% % Smooth blocks
% idx         =   s_idx(seg(1)+1:seg(2));
% mea         =   repmat(mean(X(:,idx), 2), 1, length(idx));
% C0          =   A*(X(:,idx)-mea);
% tau         =   par.tau;
% if  flag==1
%    tau      =   par.c1./(sqrt(mean(C0.^2, 2)) + e1);
%    tau      =   repmat(tau, 1, length(idx)); 
% end
% Y(:, idx)   =   A'*soft( C0, tau ) + mea;
% 
% 
% for   i  = 2:length(seg)-1
%    
%     idx  =  s_idx(seg(i)+1:seg(i+1));
%     
%     cls  =  PCA_idx(idx(1));
%     P    =   reshape(PCA_D(:, cls), b2, b2);
% 
%     mea        =   repmat(mean(X(:,idx), 2), 1, length(idx));
%     C0          =   P*(X(:,idx)-mea);
%     tau         =   par.tau;
%     if  flag==1
%         tau      =   par.c1./(sqrt(mean(C0.^2, 2)) + e1);
%         tau      =   repmat(tau, 1, length(idx)); 
%     end
%     Y(:, idx)   =   P'*soft( C0, tau ) + mea;    
% end
% 
% 
% % Output the processed image
% im_out   =  zeros(h,w);
% im_wei   =  zeros(h,w);
% k        =  0;
% for i  = 1:b
%     for j  = 1:b
%         k    =  k+1;
%         im_out(r-1+i,c-1+j)  =  im_out(r-1+i,c-1+j) + reshape( Y(k,:)', [N M]);
%         im_wei(r-1+i,c-1+j)  =  im_wei(r-1+i,c-1+j) + 1;       
%         
% %         im_out(i:h-b+i,j:w-b+j)  =  im_out(i:h-b+i,j:w-b+j) + reshape( Y(k,:)', [N M]);
% %         im_wei(i:h-b+i,j:w-b+j)  =  im_wei(i:h-b+i,j:w-b+j) + 1;
%     end
% end
% 
% im_out  =  im_out./(im_wei+eps);




% function  im_out   =  LPCA_IS_fast( im, PCA_idx, PCA_D, par, s_idx, seg, A, flag, Tau )
% 
% b        =  par.win;
% s        =  par.step;
% b2       =  b*b;
% [h  w]   =  size(im);
% 
% N       =  h-b+1;
% M       =  w-b+1;
% r       =  [1:s:N];
% r       =  [r r(end)+1:N];
% c       =  [1:s:M];
% c       =  [c c(end)+1:M];
% 
% N       =  length(r);
% M       =  length(c);
% X       =  zeros(b2, N*M);
% Y       =  X;
% 
% % For the Y component
% k    =  0;
% for i  = 1:b
%     for j  = 1:b
%         k    =  k+1;
%         blk  =  im(r-1+i,c-1+j);
%         X(k,:) =  blk(:)';
%     end
% end
% 
% 
% % Smooth blocks
% idx         =   s_idx(seg(1)+1:seg(2));
% tau         =   par.tau;
% if  flag==1
%    tau    =  Tau(:, idx); 
% end
% mea         =   repmat(mean(X(:,idx), 2), 1, length(idx));
% Y(:, idx)   =   A'*soft( A*(X(:,idx)-mea), tau ) + mea;
% 
% 
% for   i  = 2:length(seg)-1
%    
%     idx  =  s_idx(seg(i)+1:seg(i+1));
%     
%     cls  =  PCA_idx(idx(1));
%     P    =   reshape(PCA_D(:, cls), b2, b2);
% 
%     tau   =  par.tau;
%     if  flag==1
%        tau    =  Tau(:, idx); 
%     end
%     mea        =   repmat(mean(X(:,idx), 2), 1, length(idx));
%     Y(:, idx)  =   P'*soft(P*(X(:, idx)-mea), tau) + mea;
% end
% 
% 
% % Output the processed image
% im_out   =  zeros(h,w);
% im_wei   =  zeros(h,w);
% k        =  0;
% for i  = 1:b
%     for j  = 1:b
%         k    =  k+1;
%         im_out(r-1+i,c-1+j)  =  im_out(r-1+i,c-1+j) + reshape( Y(k,:)', [N M]);
%         im_wei(r-1+i,c-1+j)  =  im_wei(r-1+i,c-1+j) + 1;       
%         
% %         im_out(i:h-b+i,j:w-b+j)  =  im_out(i:h-b+i,j:w-b+j) + reshape( Y(k,:)', [N M]);
% %         im_wei(i:h-b+i,j:w-b+j)  =  im_wei(i:h-b+i,j:w-b+j) + 1;
%     end
% end
% 
% im_out  =  im_out./(im_wei+eps);

