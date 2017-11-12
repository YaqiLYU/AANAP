function S = ransac_global_similarity(data,data_orig,img1,img2)
thr_l = 0.001;
M = 500;

figure(1);
imshow([img1 img2]);
title('Ransac''s results');
hold on;
plot(data_orig(1,:),data_orig(2,:),'go','LineWidth',2);
plot(data_orig(4,:)+size(img1,2),data_orig(5,:),'go','LineWidth',2);
hold on;
pause(0.5)

for i = 1:20
    % Multi-GS内点估计
    [ ~,res,~,~ ] = multigsSampling(100,data,M,10);
    con = sum(res<=thr_l);
    [ ~, maxinx ] = max(con);
    inliers = find(res(:,maxinx)<=thr_l);
    if size(inliers) < 50
        break;
    end
    data_inliers = data(:,inliers);

    % 提取参数计算S
    x  = data_inliers(1,:); 
    y  = data_inliers(2,:); 
    x_ = data_inliers(4,:);     
    y_ = data_inliers(5,:);
    
    % 计算全局相似变换
    A = [];
    b = [];
    
    % 相似变换：平移，旋转，尺度
    for idx = 1:size(x,2)
        A = [A; x(idx) -y(idx) 1 0;
                y(idx)  x(idx) 0 1];

        b = [b;x_(idx);
               y_(idx)];
    end
    beta = A\b;
    
    S_segment{i} = [beta(1) -beta(2) beta(3);
                    beta(2)  beta(1) beta(4);
                         0        0       1];
    theta(i)     = atan(beta(2)/beta(1));
    
    % 画内点分割结果
    clr = [rand(),0,rand()];
    plot(data_orig(1,inliers),data_orig(2,inliers),...
         'o','color',clr,'LineWidth',2);
    plot(data_orig(4,inliers)+size(img1,2),data_orig(5,inliers),...
         'o','color',clr,'LineWidth',2);
    hold on;
    pause(0.5);

    % 移除内点，更新循环变量
    outliers = find(res(:,maxinx)>thr_l);
    data = data(:,outliers);
    data_orig = data_orig(:,outliers);
end

% 旋转最小的S为最优
index = find(abs(theta) == min(abs(theta)));
S = S_segment{index};
end
