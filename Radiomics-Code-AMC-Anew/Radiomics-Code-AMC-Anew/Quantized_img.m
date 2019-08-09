function [Processed_img, status, msg] = Quantized_img(img_data, mask, bins, type, alpha)
% By Bumwoo Park
% Update: 2018-01-23
% E-mail: julius0628@gmail.com


Processed_img = [];
status = 1;
msg = '';

% 0. initial checks
if sum(size(img_data) ~= size(mask))
    status = -1;
    msg = 'image and mask have different matrix size';
    return;
end

data_mask = img_data(mask);

switch type
    case 1
        % Min-Max
        min_data = min(data_mask(:));        
        max_data = max(data_mask(:));        
    case 2
        % Mean +- Std
        mean_val = mean(data_mask(:));        
        std_val = std(data_mask(:));        
        min_data = mean_val - alpha * std_val;
        max_data = mean_val + alpha * std_val;
end
data_range = min_data:(max_data-min_data)/bins:max_data;
data_range = data_range(2:end-1);

Processed_img = imquantize(img_data, data_range);
end