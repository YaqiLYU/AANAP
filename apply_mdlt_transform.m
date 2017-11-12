function [out_x, out_y] = apply_mdlt_transform(x, y, Hmdlt)

out_x = zeros(size(x));
out_y = zeros(size(y));

for i = 1:100
    for j = 1:100
        
        grididx = (i-1) * size(y,2) + j;
        T = reshape(Hmdlt(grididx,:),3,3);
        T = inv(T);
        
        in_x = x(j,i);
        in_y = y(j,i);
        
        out_x(j,i) = (T(1, 1)*in_x + T(1, 2)*in_y + T(1, 3)) ./ ...
                     (T(3, 1)*in_x + T(3, 2)*in_y + T(3, 3));

        out_y(j,i) = (T(2, 1)*in_x + T(2, 2)*in_y + T(2, 3)) ./ ...
                     (T(3, 1)*in_x + T(3, 2)*in_y + T(3, 3));
    end
end
