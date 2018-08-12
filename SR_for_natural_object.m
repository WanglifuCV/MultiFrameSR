clear; clc; close all;
addpath('BM3D');

isDisp = 1;

%% aquire data
% sequenceFile = 'data\test\LRPlatetestseq.mat';
% sequenceFile = 'data\test\LRTextseq.mat';
% sequenceFile = 'data\test\LRSuzietestseq.mat';
% sequenceFile = 'data\test\LRSuzieseq.mat';
% sequenceFile = 'data\test\LRFormanseq.mat';
sequenceFile = 'data\test\MissAmerican.mat';
[p, f, x] = fileparts(sequenceFile);
output_dir = ['result\' f '\'];
[success, message] = mkdir(output_dir);  %#ok<NASGU>

%% parameter setting
refNum = 23;   % 参考帧号
NumOfFrame = 15; % 指定总共要截取视频序列帧数

D = 10;
M = (2*D+1)^2;


% parameters
ratio = 3;  % 放大倍数
sample_org = 2; % 采用原点

iterate = 2;
window_sz = 15; % inintial window size for antialiasing filter
sigma_blur = 2;  % inintial sigma for antialiasing filter
alpha = 0.7; % decay factor for sigma of antialiasing filter
sigma_noise = 15; % initail sigma for non-local filter
beta = 0.3; % decay factor for sigma of non-local filter
weight_num = 9;
threshold_outlier = 0.2;
sigma_similar = 0.01;    % 计算块匹配相似性权重的归一化系数
sigma_similar = 2 * sigma_similar^2;
%

%% prepare for reconstruction
[LRimage refNum] = AquireSequence(sequenceFile, refNum, NumOfFrame);
[rowNum, colNum, T] = size(LRimage);
rowNumSR = rowNum * ratio; colNumSR = colNum * ratio;

Fm = zeros(2, M);
Cm = zeros(2*D+1);
k = 1;
for dy = D:-1:-D
    for dx = D:-1:-D
        Fm(:,k) = [dx; dy];
        Cm(dx + (D+1), dy + (D+1)) = k;
        k = k+1;
    end
end

im_td = zeros(rowNum, colNum, T); % save the the results of translated and decimated SR
LRimage_blur = zeros(rowNum, colNum, T); % save the antialiased observation

% Gain the initial SR result using lanczos interpolation
if isDisp
    disp('Start SR reconstruction...');
    disp('Generate the initial SR using interpolation method.');
end

Z = imresize(LRimage(:,:,refNum), ratio, 'lanczos2');
figure, imshow(Z); title('Lanczos Interpolation');
imwrite(uint8(Z*255), [output_dir 'result_Interpolated.tif'], 'tif');


Dist = zeros(rowNum, colNum, weight_num, T);
Motion = zeros(rowNum, colNum, weight_num, T);

D_M = zeros(rowNum*colNum, weight_num);
M_M = zeros(rowNum*colNum, weight_num);



%% SR reconstruction
for iter = 1 : iterate

    if isDisp
        fprintf('Start the %dth iteration...\n', iter);
    end

    % determine the sigma for blur kernel
    sigma = sigma_blur * alpha^(iter-1);
    H = fspecial('gaussian', floor(sigma*4)*2+1, sigma);

    % blur the images to reduce aliasing
    % the blured image used only for registering, not for fusion
    if isDisp
        fprintf('Anti-aliasing process... (the sigma = %f)\n', sigma);
    end
    
    if iter == 2
        sigma = 0.2;
    end
    
    if sigma > 0.3
        for i = 1:T
            LRimage_blur(:,:,i) = imfilter(LRimage(:,:,i), H, 'symmetric', 'same');
        end
    else
        LRimage_blur = LRimage;
    end

    if isDisp
        disp('Anti-aliasing Over.')
    end

    % Translate ,Decimate and Blur the SR of refrence frame
    if isDisp
        disp('Gernerate the translated and decimated examples of SR...');
    end

    k = 1;
    for dy = D:-1:-D
        for dx = D:-1:-D
            [xx, yy] = meshgrid(1:size(Z,2), 1:size(Z,1));
            xx2 = xx - dx; yy2 = yy - dy;
            im_shift = interp2(xx,yy, Z, xx2, yy2, 'cubic', NaN);
            I = isnan(im_shift);
            im_shift(I) = Z(I);
            im_td(:,:,k) = im_shift(sample_org:ratio:end, sample_org:ratio:end);
            if sigma > 0.3
                im_td(:,:,k) = imfilter(im_td(:,:,k), H, 'symmetric', 'same');
            end

            if isDisp
                fprintf('%d...', k);
            end

            k = k+1;

        end
    end

    if isDisp
        fprintf('\nGernerate the translated and decimated examples of SR Over.\n')
    end


    % determine the size of matching window for registering
    w_sz = window_sz - 4*(iter-1);
    if w_sz < 3
        w_sz = 3;
    end
    
    if iter == 2
        w_sz = 3;
        threshold_outlier = 1;
    end

    % compute threshold
    H = ones(w_sz);
    search_range = ceil((ratio-1)/2);  % 选择移动较小的样本计算，确定范围。
    Thres = zeros(rowNum*colNum, 1);
    for m = 1:M
        if Fm(2, m) ~= 0 || Fm(1, m) ~= 0
            if abs(Fm(2, m))<=search_range && abs(Fm(1, m))<=search_range
                Err = (im_td(:,:,m) - LRimage_blur(:,:,refNum)).^2;
                Err = conv2(Err, H, 'same');
                Err = Err(:);
                indx = Thres < Err;
                Thres(indx) = Err(indx);
            end
        end
    end

    Thres = Thres * threshold_outlier;
    Thres = repmat(Thres, [1 weight_num]);

    for t = 1:T
        D_M(:,:) = Inf;
        M_M(:,:) = 0;
        if t ~= refNum
            for m = 1:M
                Err = (im_td(:,:,m) - LRimage_blur(:,:,t)).^2;
                Err = conv2(Err, H, 'same');
                Err = Err(:);

                for i=1:weight_num
                    indx = D_M(:,i) > Err;
                    if sum(indx)~=0
                        for j = weight_num : -1: i+1
                            D_M(indx, j) = D_M(indx, j-1);
                            M_M(indx, j) = M_M(indx, j-1);
                        end
                        D_M(indx, i) = Err(indx);
                        M_M(indx, i) = m;
                        Err(indx) = Inf;
                    end
                end
            end

            % 去除 outlier
            C = D_M > Thres;
            D_M(C) = Inf;
            M_M(C) = 0;

        else
            D_M(:,1) = 0;
            D_M(:,2:weight_num) = Inf;
            M_M(:,1) = Cm(D+1,D+1);
        end

        Dist(:,:,:,t) = reshape(D_M, rowNum, colNum, weight_num);
        Motion(:,:,:,t) = reshape(M_M, rowNum, colNum, weight_num);

    end


    %fusion
    val = zeros(rowNumSR, colNumSR);
    normal = zeros(rowNumSR, colNumSR);
    for t = 1:T
        for k = 1:colNum
            for l = 1:rowNum
                for  i = 1:weight_num
                    w = exp(-Dist(l, k, i, t)/(sigma_similar * w_sz^2));
                    if w > 0
                        m = Motion(l, k, i, t);
                        x = (k-1)*ratio + sample_org - Fm(1, m);
                        y = (l-1)*ratio + sample_org - Fm(2, m);
                        if (x>=1 && x<=colNumSR && y>=1 && y<=rowNumSR )
                            val(y, x) = val(y, x) + w * LRimage(l,k,t);
                            normal(y, x) = normal(y, x) + w;
                        end
                    end
                end
            end
        end
    end

    C = normal > 0;
    Z(C) = val(C)./normal(C);

    figure; imshow(Z); title(['the ' num2str(iter) ' iteration SR result']);
    imwrite(uint8(Z*255), [output_dir 'result_SR_' num2str(iter) '.tif'], 'tif');
    
    % Non-Local Filter
%     sigma_n = sigma_noise * beta^(iter - 1);
%       sigma_n = 25;
    sigma_n = sigma_noise - 2* iter;
    if isDisp
        fprintf('Non-local filtering (the sigma = %f)...\n', sigma_n);
    end

    [NA, ZI] = BM3D(1, Z, sigma_n);
    
    % modify the result by block sparsity
    Zs = (ZI-Z).^2;
    ths = max(Zs(:)) * 0.02*((iterate - iter+1)*2);
    Indx = Zs > ths;
    Z(Indx) = ZI(Indx); 
    Z(sample_org:ratio:end, sample_org:ratio:end) = LRimage(:,:,refNum);
    
    if isDisp
        disp('Non-local filtering Over');
    end
    
    figure; imshow(Z); title(['the ' num2str(iter) ' iteration filterd SR result']);
    imwrite(uint8(Z*255), [output_dir 'result_filtered_SR_' num2str(iter) '.tif'], 'tif');
    
    if isDisp
        fprintf('the %dth iteration Over.\n', iter);
    end

end



if isDisp
    disp('SR reconstruction Over.');
end




