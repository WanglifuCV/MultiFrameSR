clear; clc; close all;
addpath('BM3D');

isDisp = 1;

% aquire data
sequenceFile = 'data\test\LRPlatetestseq.mat';
% sequenceFile = 'data\test\LRTextseq.mat';
% sequenceFile = 'data\test\LRSuzietestseq.mat';
% sequenceFile = 'data\test\LRSuzieseq.mat';
% sequenceFile = 'data\test\LRFormanseq.mat';
[p, f, x] = fileparts(sequenceFile);
output_dir = ['result\' f '\']; 
[success, message] = mkdir(output_dir);  %#ok<NASGU>

%% parameter setting
refNum = 5;   % 参考帧号
NumOfFrame = 9; % 指定总共要截取视频序列帧数

% parameters
ratio = 3;  % 放大倍数
sample_org = 2; % 采用原点

D = 1;
M = (2*D+1)^2;

iterate = 2;
window_sz = 15; % inintial window size for antialiasing filter
sigma_blur = 1.35;  % inintial sigma for antialiasing filter
alpha = 0.5; % decay factor for sigma of antialiasing filter
sigma_noise = 25; % initail sigma for non-local filter
beta = 0.3; % decay factor for sigma of non-local filter
sigma_similar = 0.01;    % 计算块匹配相似性权重的归一化系数
sigma_similar = 2 * sigma_similar^2;


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

%% SR reconstruction
for iter = 1 : iterate
    
    if isDisp
        fprintf('Start the %dth iteration...\n', iter);
    end

    val = zeros(rowNumSR, colNumSR);
    normal = zeros(rowNumSR, colNumSR);
    
    % determine the sigma for blur kernel
    sigma = sigma_blur * alpha^(iter-1);
    H = fspecial('gaussian', floor(sigma*4)*2+1, sigma);

    % blur the images to reduce aliasing
    % the blured image used only for registering, not for fusion
    if isDisp
        fprintf('Anti-aliasing process... (the sigma = %f)\n', sigma);
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
            im_td(:,:,k) = im_shift(2:3:end, 2:3:end);
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

    % detect points fit for registering
    Map_corner = HarrisDetector(LRimage_blur(:,:,refNum), w_sz);
    
    if isDisp
        disp('Registering and Fusion process...');
    end
    
    Motion = []; Probability = [];
    for t = 1:T
        % Register
        if isDisp
            fprintf('Registering for the %d frame...\n', t);
        end
        
        if t ~= refNum
            [Motion Probability] = RegisterForAFrame(LRimage_blur(:,:,t), im_td, Map_corner, sigma_similar, w_sz, D, Fm, Cm);
        else
            Motion = zeros(rowNum, colNum, 2);
            Probability = ones(rowNum, colNum);
        end
        
        if isDisp
            fprintf('Registering for the %d frame Over.\nFusing for the %d frame...\n', t, t);
        end
        
        % Fusion
        for k = 1:colNum
            for l = 1:rowNum
                w = Probability(l, k);
                if w > 0
%                 if w > 10^(-4)
                    x = (k-1)*ratio + sample_org - Motion(l, k, 1);
                    y = (l-1)*ratio + sample_org - Motion(l, k, 2);

                    if (x>=1 && x<=colNumSR && y>=1 && y<=rowNumSR )
                        val(y, x) = val(y, x) + w * LRimage(l, k, t);
                        normal(y, x) = normal(y, x) + w;
                    end
                end
            end
        end
        if isDisp
            fprintf('Fusing for the %d frame Over.\n', t);
        end
    end

    Idx = normal > 0;
    Z(Idx) = val(Idx)./normal(Idx);
    
    figure; imshow(Z); title(['the ' num2str(iter) ' iteration SR result']);
    imwrite(uint8(Z*255), [output_dir 'result_SR_' num2str(iter) '.tif'], 'tif');
    
    % Non-Local Filter
    sigma_n = sigma_noise * beta^(iter - 1);
    
    if isDisp
        fprintf('Non-local filtering (the sigma = %f)...\n', sigma_n);
    end
    
    [NA, Z] = BM3D(1, Z, sigma_n);
%     Z = medfilt2(Z, [3 3]);

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
