function [ H_output ] = homography_linearization( H, point, anchor_points )
vega = 5;

% Obtain kernel: Student’s t-weighting 
alpha = (1 + pdist2(point,anchor_points)./vega).^(-(vega+1)/2);
alpha = alpha./sum(alpha);

% Linearization using Taylor series
H_output = zeros(3,3);

for i = 1:size(anchor_points,1)  
    % 依次在每个锚点展开
    A = taylor_series( H, anchor_points(i,:));
    
    % 加权平均
    H_output = H_output + alpha(i).*A;
end

end

function [ A ] = taylor_series( H, p)
% The first two terms of Taylor series provide the best linearization of H
% H:3*3    p:1*2 

h1 = H(1,1); h2 = H(1,2); h3 = H(1,3);
h4 = H(2,1); h5 = H(2,2); h6 = H(2,3);
h7 = H(3,1); h8 = H(3,2); h9 = H(3,3);
x1 = p(1,1);x2 = p(1,2);

dy1dx1 = h1/(h9 + h7*x1 + h8*x2) - (h7*(h3 + h1*x1 + h2*x2))/(h9 + h7*x1 + h8*x2)^2;
dy1dx2 = h2/(h9 + h7*x1 + h8*x2) - (h8*(h3 + h1*x1 + h2*x2))/(h9 + h7*x1 + h8*x2)^2;
dy2dx1 = h4/(h9 + h7*x1 + h8*x2) - (h7*(h6 + h4*x1 + h5*x2))/(h9 + h7*x1 + h8*x2)^2;
dy2dx2 = h5/(h9 + h7*x1 + h8*x2) - (h8*(h6 + h4*x1 + h5*x2))/(h9 + h7*x1 + h8*x2)^2;

y1 = (h3 + h1*x1 + h2*x2)/(h9 + h7*x1 + h8*x2);
y2 = (h6 + h4*x1 + h5*x2)/(h9 + h7*x1 + h8*x2);

A = [dy1dx1 dy1dx2 (y1 - x1*dy1dx1 - x2*dy1dx2);
     dy2dx1 dy2dx2 (y2 - x1*dy2dx1 - x2*dy2dx2);
          0      0                           1];
end
