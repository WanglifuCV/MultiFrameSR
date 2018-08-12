%-------------------------------
% PCA dictionaries and AR model
% training
% Modified, 09/May/2010, Weisheng Dong
% wsdong@mail.xidian.edu.cn
%-------------------------------
% clc;
clear;

TD_path      =   'Data/Data_set/TD2';

% Set the parameters
cls_num     =   200;
b           =   7;
ws          =   3;
delta       =   4.4;
psf         =   fspecial('gauss', 7, 1.5);
num_blks    =   400*(b^2)*40;

fpath       =   fullfile(TD_path, '*.tif');
im_dir      =   dir(fpath);
im_num      =   length(im_dir);

X     =  zeros(0);
Y     =  zeros(0);
for  i  =  1 : im_num

    im         =   double( imread(fullfile(TD_path, im_dir(i).name)) );
    [Py Px]    =   Get_patches( im, b, psf );
    X          =   [X, Px];
    Y          =   [Y, Py];
end

% Select the smooth blocks
m     =  mean(Y);
d     =  (Y-repmat(m, size(Y,1), 1)).^2;
v     =  sqrt( mean( d ) );
[a, idx]  =  find( v>=delta );
[a, idx2] =  find( v<delta );

X0   =   X(:, idx2);
[P0, mx]   =  getpca(X0);
clear X0 Px Py a d m v idx2;

X    =   X(:, idx);
Y    =   Y(:, idx);
N    =   size(X, 2);
P    =   randperm(N)';
num  =   min(N, num_blks);
P    =   P(1:num);
X    =   X(:, P);
Y    =   Y(:, P);
clear P;

fprintf('\n\nLearning On Training set 1...\n');
[centroids, cls_idx, s_idx, seg, cls_num]   =  Clustering( Y, cls_num );


%Training PCA dictionaries and AR models
PCA_D    =  zeros(b^4, cls_num);
AR_D     =  zeros(ws*ws-1, cls_num);
for  i  =  1 : length(seg)-1
   
    idx    =   s_idx(seg(i)+1:seg(i+1));    
    cls    =   cls_idx(idx(1));    
    X1     =   X(:, idx);
    
    % Compute the PCA basis
    [P, mx]   =  getpca(X1);    
    PCA_D(:,cls)    =  P(:);
    
    % Compute the AR models
    AR    =  Get_AR( X1, ws );
    AR_D(:, cls)    =  AR;
end

PCA_D    =  [P0(:) PCA_D];
save Data/Lib/PCA_AR_TD2_c200 centroids PCA_D AR_D;


% %==========================================================================
% %Training PCA dictionaries and AR models
% Var_D    =  zeros(b^2, cls_num);
% 
% for  i  =  1 : length(seg)-1
%    
%     idx    =   s_idx(seg(i)+1:seg(i+1));    
%     cls    =   cls_idx(idx(1));    
%     X1     =   X(:, idx);
%     
%     % Compute the PCA basis
%     [P, mx]   =  getpca(X1);    
%     var       =  P*(X1-mx(:, ones(1,length(idx))));
%     var       =  sqrt(mean(var.^2, 2));
%     Var_D(:,cls)   =  var;
%     
% end
% 
% save Data/Lib/PCA_AR_TD2 centroids PCA_D AR_D Var_D;



