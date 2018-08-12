function  [cent, cls_idx, s_idx, seg, cls_num]   =  Enhanced_clustering( X, cls_num )

b2      =  size(X, 1);
L       =  size(X, 2);
itn     =  12;

P      =  randperm(L);
P2     =  P(1:cls_num);
vec    =  X(:,P2(1:end));


for i = 1 : itn
    
    D         =  zeros(1, cls_num);
    dis       =  ones(1, L, 'single')*10000000000;
    cls_idx   =  zeros(L,1);
    
    for  j  = 1:L
       
       cb   =  repmat( X(:,j), 1, cls_num );
       d    =  sum( (cb-vec).^2 );
       [v idx]  =  min(d);
       dis(j)   =  v;
       cls_idx(j)  =  idx;
    end
    
    [s_idx, seg]   =  Proc_cls_idx( cls_idx );
    
    for  j  = 1:length(seg)-1   
        idx  =  s_idx(seg(j)+1:seg(j+1));
        cls  =  cls_idx(idx(1));        
        D(cls)   =  sum( dis(idx) );
    end
    mD    =  mean( D );
    U     =  D./mD;
    Sp    =  zeros(cls_num,1);
    Un    =  ones(1,cls_num);
    
    % Shifting the centroids
    for  k  = 1 : cls_num
        
        if  U(k)<1 && Sp(k)==0
            idx_k  =  find(cls_idx==k);
            
            if isempty( idx_k )
                l  =  0;
            else
                % find similar vector of kth centroid
                v2     =  vec;
                v2(:,k)=  [];
                cb     =  repmat(vec(:,k), 1, cls_num-1);
                d      =  sum( (cb-v2).^2 );
                [d l]  =  min(d);
                
                if Sp(l)==1
                    l  =  -1;
                else
                    
                    idx_l  =  find( cls_idx==l );
                end                
            end;
            
            if l<0
                continue;
            end
            
            % Calculate new 
            idx_kl  =  [idx_l; idx_k];
            vl      =  mean( X(:, idx_kl), 2 );
            dl      =  sum(sum( (X(:,idx_kl) - repmat(vl, 1, size(idx_kl,1))).^2 ));
            
            
            % select a cell whose Ui>1
            s1    =  find((U.*Un)>1);
            P     =  U(s1)./sum(U(s1));
            P     =  cumsum(P);            
            t     =  find(rand<=P);
            
            if isempty(t)
                break;
            end
            
            p  =  t(1);                  
            idx_p  =  find( cls_idx==p );
            
            % Generate two centroids
            v_p    =  X(:, idx_p);
            min_p  =  min( v_p,[],2 );
            max_p  =  max( v_p,[],2 );
            vk     =  min_p + (max_p - min_p)./4;
            vp     =  min_p + (max_p - min_p).*0.75;
            
            
            Y    =  X(:,idx_p);
            L2   =  size(Y,2);
            dd   =  zeros(2, L2);
            for  j  =  1 : 8
                
                ck   =  repmat( vk, 1, L2 );
                cp   =  repmat( vp, 1, L2 );
                dd(1,:)   =  sum( (Y-ck).^2 );
                dd(2,:)   =  sum( (Y-cp).^2 );
                [vd it]   =  min( dd );
                
                vk   =  mean(Y(:, it==1),2);
                vp   =  mean(Y(:, it==2),2);                
            end
            
            if l==0
                d_old    =  D(p);
                d_new    =  sum(vd);
            else
                d_old    =  D(l) + D(k) + D(p);
                d_new    =  dl + sum(vd);
            end
                        
            if  d_new < d_old
               cls_idx( idx_kl )   =  l;
               
               ik   =  find(it==1);
               ip   =  find(it==2);
               cls_idx( idx_p(ik) )   =  k;
               cls_idx( idx_p(ip) )   =  p;
               
               if l>0
                   D(l)   =  dl;
                   U(l)   =  D(l)/mD;
               end
               D(k)   =  sum(vd( ik ));
               D(p)   =  sum(vd( ip ));
               U(k)   =  D(k)/mD;
               U(p)   =  D(p)/mD;
               
               Un(l)  =  0;
               Sp(p)  =  1;            
            end
        end
    end
    [s_idx, seg]   =  Proc_cls_idx( cls_idx );
    
    cls_cnt    =  zeros(cls_num, 1);
    
    for  j  = 1:length(seg)-1   
         idx  =  s_idx(seg(j)+1:seg(j+1));
         cls  =  cls_idx(idx(1));
         vec(:,cls)  =  mean( X(:,idx), 2 );
         
         cls_cnt( cls )  =  length(idx);
    end
    
    mse = sum(D)/L/b2;
    disp( sprintf('clustering %d th loop, mse = %f', i, mse) );
end

save CLS_CNT  cls_cnt;
[val ind]  =  min( cls_cnt );       % Remove these classes with little samples

% remove the clusters
while val<300
    vec(:,ind)    =  [];
    cls_num       =  cls_num - 1;
    cls_cnt(ind)  =  [];
    [val  ind]    =  min(cls_cnt);
end        

cls_idx   =  zeros(L,1);
for  j  = 1:L
       
   cb   =  repmat( X(:,j), 1, cls_num );
   d    =  sum( (cb-vec).^2 );
   [v idx]  =  min(d);
   cls_idx(j)  =  idx;
end

[s_idx, seg]   =  Proc_cls_idx( cls_idx );
    
cls_cnt    =  zeros(cls_num, 1);
    
for  j  = 1:length(seg)-1   
     idx  =  s_idx(seg(j)+1:seg(j+1));
     cls  =  cls_idx(idx(1));
     vec(:,cls)  =  mean( X(:,idx), 2 );
         
     cls_cnt( cls )  =  length(idx);
end

cls_num
cent   =  vec;


