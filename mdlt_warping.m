function warped_img = mdlt_warping(height, width, img, Hmdlt, offset, X, Y)

warped_img = zeros(height,width,3);
[imgh,imgw,~] = size(img);

for xidx= 1:width
    for yidx = 1:height
        
        x_grididx = min(find(X >= xidx,1));
        y_grididx = min(find(Y >= yidx,1));
        
        grididx = (x_grididx-1) * size(Y,2) + y_grididx;
        h = reshape(Hmdlt(grididx,:),3,3);
        
        x = xidx - offset(1) + 1;
        y = yidx - offset(2) + 1;
        
        % [_x;_y;1] = H * [x;y;1]   
        posz =  h(3,1) * x + h(3,2) * y + h(3,3) ;
        posx = (h(1,1) * x + h(1,2) * y + h(1,3)) ./ posz;
        posy = (h(2,1) * x + h(2,2) * y + h(2,3)) ./ posz;

        if (posx>=1)&&(posx<=imgw)&&(posy>=1)&&(posy<=imgh)
%             warped_img(yidx,xidx,:) = img(posy,posx,:);
            
            pix = [posy posx];
            float_Y=pix(1)-floor(pix(1)); 
            float_X=pix(2)-floor(pix(2));
            
            pix_up_left=[floor(pix(1)) floor(pix(2))];
            pix_up_right=[floor(pix(1)) ceil(pix(2))];
            pix_down_left=[ceil(pix(1)) floor(pix(2))];
            pix_down_right=[ceil(pix(1)) ceil(pix(2))];
            
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
