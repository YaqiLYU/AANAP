function [output_canvas] = image_blending_linear(warped_img1,warped_img2)

    c1omask{1} = imfill(im2bw(uint8(warped_img1), 0),'holes');
    c1omask{2} = imfill(im2bw(uint8(warped_img2), 0),'holes');

    c1out{1} = double(warped_img1);
    c1out{2} = double(warped_img2);
    
    for i = 1 : 2
        if i == 1
            out = c1out{1};
            out_mask = c1omask{1};
            % center of out
            [r, c] = find(c1omask{1}(:, :, 1));
            out_center = [mean(r) mean(c)];
        else % blend out and c1out{i}
            % center of c1out{i}
            [r, c] = find(c1omask{i}(:, :, 1));
            out_i_center = [mean(r) mean(c)];
            % compute weighted mask
            vec = out_i_center - out_center; % vector from out_center to out_i_center
            intsct_mask = c1omask{i}(:, :, 1) & out_mask(:, :, 1); % 1 channel

            [r, c] = find(intsct_mask(:, :, 1));
            idx = sub2ind(size(c1omask{i}(:, :, 1)), r, c);
            out_wmask = zeros(size(c1omask{i}(:, :, 1)));
            proj_val = (r - out_center(1))*vec(1) + (c- out_center(2))*vec(2); % inner product
            out_wmask(idx) = (proj_val - (min(proj_val)+(1e-3))) / ...
                             ((max(proj_val)-(1e-3)) - (min(proj_val)+(1e-3))); % weight map (of overlapped area) for c1out{i}, 1 channel
            % blending
            mask1 = out_mask(:, :, 1)&(out_wmask==0);
            mask2 = out_wmask;
            mask3 = c1omask{i}(:, :, 1)&(out_wmask==0);
            mask1 = cat(3, mask1, mask1, mask1); mask2 = cat(3, mask2, mask2, mask2); mask3 = cat(3, mask3, mask3, mask3);
            out = out.*(mask1+(1-mask2).*(mask2~=0)) + c1out{i}.*(mask2+mask3);
            % update
            out_mask = out_mask | c1omask{i};
            out_center = out_i_center; % update out_center by assign center of c1out{i}
        end
    end
    
    output_canvas = uint8(out);
end
