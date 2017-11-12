%**************************************************%
% An implementation of AANAP below:
%
% Lin C C, Pankanti S U, Ramamurthy K N, et al.
% Adaptive as-natural-as-possible image stitching [C]
% 2015 IEEE Conference on Computer Vision and Pattern Recognition (CVPR). 
% IEEE, 2015: 1155-1163.
% 
% Just for fun! For commercial purposes, please contact the author.
% by YaqiLYU
% https://www.zhihu.com/question/34535199/answer/135169187
%**************************************************%

close all

run('vlfeat-0.9.20/toolbox/vl_setup');
addpath('modelspecific');
addpath('mexfiles');
addpath('multigs');

%% Global options
% 0 - Bilinear interpolation, implementation by MATLAB，slower but better
% 1 - Nearest neighbor interpolation,implementation by C++, Faster but worse
fast_stitch = 1;    
img_n = 2;  % only support two image stitching
in_name = cell(img_n,1);
in_name{1} = 'images\railtracks_01.jpg';
in_name{2} = 'images\railtracks_02.jpg';
img_n = size(in_name, 1);

gamma = 0;
sigma = 12.5;

%% load and preprocessing
I = cell(img_n, 1);
for i = 1 : img_n
    I{i} = imread(in_name{i});
end

max_size = 1000 * 1000;
imgw = zeros(img_n, 1);
imgh = zeros(img_n, 1);

for i = 1 : img_n
    if numel(I{i}(:, :, 1)) > max_size
        I{i} = imresize(I{i}, sqrt(max_size / numel(I{i}(:, :, 1))));
    end

    imgw(i) = size(I{i}, 2);
    imgh(i) = size(I{i}, 1);
end

img1 = I{1};
img2 = I{2};
img2 = imresize(img2,size(img1,1)/size(img2,1));

figure(4),
imshow(img1,[]);
hold on;
pause(0.3);
imshow(img2,[]);
pause(0.3);

%% User defined parameters for APAP
clear global;
global fitfn resfn degenfn psize numpar
fitfn = 'homography_fit';
resfn = 'homography_res';
degenfn = 'homography_degen';
psize   = 4;
numpar  = 9;

M     = 500;
thr_g = 0.1;

if fast_stitch
    C1    = 100;
    C2    = 100;
else
    C1    = 200;
    C2    = 200;
end

%% SIFT keypoint detection and matching.
[ kp1,ds1 ] = vl_sift(single(rgb2gray(img1)),'PeakThresh', 0,'edgethresh',500);
[ kp2,ds2 ] = vl_sift(single(rgb2gray(img2)),'PeakThresh', 0,'edgethresh',500);
[match_idxs, scores] = vl_ubcmatch(ds1,ds2);
f1 = kp1(:,match_idxs(1,:));
f2 = kp2(:,match_idxs(2,:));

%% Normalise point distribution and Outlier removal with Multi-GS RANSAC.
% (x1;y1;1;x2;y2;1)
data_orig = [ kp1(1:2,match_idxs(1,:)) ; ones(1,size(match_idxs,2)) ;
              kp2(1:2,match_idxs(2,:)) ; ones(1,size(match_idxs,2)) ];
[ dat_norm_img1,T1 ] = normalise2dpts(data_orig(1:3,:));
[ dat_norm_img2,T2 ] = normalise2dpts(data_orig(4:6,:));
data_norm = [ dat_norm_img1 ; dat_norm_img2 ];

% Multi-GS
% rng(0);
[ ~,res,~,~ ] = multigsSampling(100,data_norm,M,10);
con = sum(res<=thr_g);
[ ~, maxinx ] = max(con);
inliers = find(res(:,maxinx)<=thr_g);

%% Global homography (H) again.
[ Hl,A,D1,D2 ] = feval(fitfn,data_norm(:,inliers));% 输入评估输出
Hg = T2\(reshape(Hl,3,3)*T1);% H去归一化
Hg = Hg / Hg(3,3)

%% Compute Global similarity
S = ransac_global_similarity(data_norm(:,inliers),data_orig(:,inliers),img1,img2);
S = T2\(S*T1)

%% Obtaining size of canvas (using global Homography).
TL = Hg\[1;1;1];
TL = round([ TL(1)/TL(3) ; TL(2)/TL(3) ]);
BL = Hg\[1;size(img2,1);1];
BL = round([ BL(1)/BL(3) ; BL(2)/BL(3) ]);
TR = Hg\[size(img2,2);1;1];
TR = round([ TR(1)/TR(3) ; TR(2)/TR(3) ]);
BR = Hg\[size(img2,2);size(img2,1);1];
BR = round([ BR(1)/BR(3) ; BR(2)/BR(3) ]);

% Canvas size.
cw = max([1 size(img1,2) TL(1) BL(1) TR(1) BR(1)]) - min([1 size(img1,2) TL(1) BL(1) TR(1) BR(1)]) + 1;
ch = max([1 size(img1,1) TL(2) BL(2) TR(2) BR(2)]) - min([1 size(img1,1) TL(2) BL(2) TR(2) BR(2)]) + 1;

% Offset for left image.
off = [ 1 - min([1 size(img1,2) TL(1) BL(1) TR(1) BR(1)]) + 1 ;
        1 - min([1 size(img1,1) TL(2) BL(2) TR(2) BR(2)]) + 1 ];

%% Generate anchor points in the boundary,20 in each size, 80 in total
anchor_points = [];
anchor_num = 20;
hx = linspace(1,size(img1,2),anchor_num);
hy = linspace(1,size(img1,1),anchor_num);

for i = 1:anchor_num
    anchor_points = [anchor_points;
                     0, round(hy(i))];    
    anchor_points = [anchor_points;
                     size(img1,2), round(hy(i))];   
    anchor_points = [anchor_points;
                     round(hx(i)), 0];    
    anchor_points = [anchor_points;
                     round(hx(i)), size(img1,1)];
end

%% Image stitching with global homography (H).
warped_img1 = uint8(zeros(ch,cw,3));
warped_img1(off(2):(off(2)+size(img1,1)-1),off(1):(off(1)+size(img1,2)-1),:) = img1;
warped_img2 = imagewarping(double(ch),double(cw),double(img2),Hg,double(off));
warped_img2 = reshape(uint8(warped_img2),size(warped_img2,1),size(warped_img2,2)/3,3);
linear_hom = imageblending(warped_img1,warped_img2);

figure(2),imshow(linear_hom);
title('Image Stitching with global homography');
hold on;
pause(0.3);
plot(anchor_points(:,1),anchor_points(:,2)+ off(2),'go','LineWidth',2);
pause(0.3);

%% Compute weight for Integration
% (x,y)：K_min -> K_1 -> K_2 -> K_max
Or = [size(img1,2)/2;size(img1,1)/2];
Ot = Hg\[size(img2,2)/2;size(img2,1)/2;1];
Ot = [Ot(1)/Ot(3);Ot(2)/Ot(3)];

% sovel linear problem
k = (Ot(2) - Or(2))/(Ot(1) - Or(1));
b = Or(2) - k * Or(1);

K_min(1) = min([TL(1) BL(1) TR(1) BR(1)]);
K_max(1) = max([TL(1) BL(1) TR(1) BR(1)]);
K_1(1) = size(img1,2);
K_2(1) = K_1(1) + (K_max(1) - K_1(1))/2;

K_min(2) = k * K_min(1) + b;
K_max(2) = k * K_max(1) + b;
K_1(2) = k * K_1(1) + b;
K_2(2) = k * K_2(1) + b;

% Image keypoints coordinates
% 特征点的真实坐标在图1中
Kp = [data_orig(1,inliers)' data_orig(2,inliers)'];

[ X,Y ] = meshgrid(linspace(1,cw,C1),linspace(1,ch,C2));

% Mesh (cells) vertices' coordinates.
% 全景图的真实坐标
Mv = [X(:)-off(1), Y(:)-off(2)];

% Perform Moving DLT
fprintf('  Moving DLT main loop...');tic;
Ht = zeros(size(Mv,1),9);
Hr = zeros(size(Mv,1),9);

%% Moving DLT and Similarly Extrapolate (projective).
% 重叠区域：局部单应性和全局相似变换平滑过渡
% 非重叠区域：线性化局部单应性(仿射)和全局相似变换平滑过渡

for i = 1:size(Mv,1)
    % Obtain kernel: Gaussian weighting
    Gki = exp(-pdist2(Mv(i,:),Kp)./sigma^2);   
    Wi = max(gamma,Gki); 

    v = wsvd(Wi,A);
    Hl = reshape(v,3,3)'; 
    
    % De-condition and De-normalize
    Hl = D2\Hl*D1;
    Hl = T2\Hl*T1;
    Hl = Hl / Hl(3,3);
    
    P = [Mv(i,1);Mv(i,2)];
    [a, b] = compute_weight(P, K_min, K_max);
    [c, d] = compute_weight(P, K_1, K_2);
    
    % Homography linearization in the non-overlapping region
    if Mv(i,1) >= size(img1,2)
        Hl_linear = homography_linearization(Hl, Mv(i,:), anchor_points);
        Hl = c.*Hl_linear + d.*Hl;
    end

    % Interpolate smoothly and similarly extrapolate
    h_target = a.*S + b.*Hl;
    Ht(i,:) = h_target(:);
    
    h_refer = Hl\h_target; % h_target*inv(Hl) in paper, ERROR!
    if Mv(i,1) > size(img1,2)
        h_refer = zeros(3,3);
    end
    Hr(i,:) = h_refer(:);
end
fprintf('done (%fs)\n',toc);

%% Warping images with Moving DLT
if fast_stitch
    [warped_img1] = imagewarping(double(ch),double(cw),double(img1),Hr,double(off),X(1,:),Y(:,1)');
    warped_img1 = reshape(uint8(warped_img1),size(warped_img1,1),size(warped_img1,2)/3,3);
    figure(3),imshow(uint8(warped_img1));
    hold on;
    pause(0.5)
    
    [warped_img2] = imagewarping(double(ch),double(cw),double(img2),Ht,double(off),X(1,:),Y(:,1)');
    warped_img2 = reshape(uint8(warped_img2),size(warped_img2,1),size(warped_img2,2)/3,3);
    imshow(uint8(warped_img2));
    hold on;
    pause(0.3)
    
    pano = image_blending_average(warped_img1,warped_img2);
    imshow(pano);
    title('Adaptive As-Natural-As-Possible Image Stitching');
else
    warped_img1 = mdlt_warping(ch,cw,img1,Hr,off,X(1,:),Y(:,1)');
    figure(3),imshow(uint8(warped_img1));
    hold on;
    pause(0.3)
    
    warped_img2 = mdlt_warping(ch,cw,img2,Ht,off,X(1,:),Y(:,1)');
    imshow(uint8(warped_img2));
    hold on;
    pause(0.3)
    
    pano = image_blending_linear(warped_img1,warped_img2, [off(2)+1;off(2)+size(img1,1)-1;size(img1,2)]);
    imshow(pano);
    title('Adaptive As-Natural-As-Possible Image Stitching');
end
