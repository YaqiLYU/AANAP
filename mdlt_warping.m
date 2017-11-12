%-------------------------------------------------------------------------
% 局部单应性的图像绘制函数
% width * height：目标图像的宽和高(由全局单应性方法估计得到)
% img：输入图像，即源图像  warped_img：输出图像，即目标图像
% Hmdlt：局部单应性矩阵列表，cell数量 * 9
% offset: 目标图像在拼接后图像中的偏移量，防止绘制后图像超出边界
% X, Y：网格cell的宽高位置，已经过offset偏移
% 保持一致，x与X都表示宽width,y和Y都表示高height
% 网格划分Mv与网格对应的局部单应性矩阵格式为[width,height]
% by YaqiLYU
%-------------------------------------------------------------------------
function warped_img = mdlt_warping(height, width, img, Hmdlt, offset, X, Y)

warped_img = zeros(height,width,3);
[imgh,imgw,~] = size(img);

% 遍历目标图像中的每个像素点(逆投影，可以避免空洞和重叠问题)Matlab列优先
for xidx= 1:width
    for yidx = 1:height
        
        % 当前像素点坐标的网格位置
        x_grididx = min(find(X >= xidx,1));
        y_grididx = min(find(Y >= yidx,1));
        
        % 取出对应网格的局部单应性矩阵
        grididx = (x_grididx-1) * size(Y,2) + y_grididx;
        h = reshape(Hmdlt(grididx,:),3,3);
        
        % 补偿网格平移变换
        x = xidx - offset(1) + 1;
        y = yidx - offset(2) + 1;
        
        % 目标坐标(xidx,yidx)投影到源图像坐标(posx,posy)
        % [_x;_y;1] = H * [x;y;1]   
        posz =  h(3,1) * x + h(3,2) * y + h(3,3) ;
        posx = (h(1,1) * x + h(1,2) * y + h(1,3)) ./ posz;
        posy = (h(2,1) * x + h(2,2) * y + h(2,3)) ./ posz;

        % 判断绘制后坐标是否落在源图像的画布中
        if (posx>=1)&&(posx<=imgw)&&(posy>=1)&&(posy<=imgh)
%             warped_img(yidx,xidx,:) = img(posy,posx,:);
            
            pix = [posy posx];
            float_Y=pix(1)-floor(pix(1)); 
            float_X=pix(2)-floor(pix(2));
            
            % 四个相邻的点
            pix_up_left=[floor(pix(1)) floor(pix(2))];
            pix_up_right=[floor(pix(1)) ceil(pix(2))];
            pix_down_left=[ceil(pix(1)) floor(pix(2))];
            pix_down_right=[ceil(pix(1)) ceil(pix(2))];
            
            %计算临近四个点的权重
            value_up_left   = (1-float_X) * (1-float_Y);
            value_up_right  =     float_X * (1-float_Y);
            value_down_left = (1-float_X) * float_Y;
            value_down_right=     float_X * float_Y;
                                                            
            warped_img(yidx,xidx,:) = ...
               value_up_left*img(   pix_up_left(1),   pix_up_left(2),:)+ ...
              value_up_right*img(  pix_up_right(1),  pix_up_right(2),:)+ ...
             value_down_left*img( pix_down_left(1), pix_down_left(2),:)+ ...
            value_down_right*img(pix_down_right(1),pix_down_right(2),:);

        else
            warped_img(yidx,xidx,:) = [0,0,0];
        end
    end
end

end
