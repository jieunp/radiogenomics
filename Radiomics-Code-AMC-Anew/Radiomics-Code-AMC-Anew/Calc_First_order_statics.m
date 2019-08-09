function [stat, status, msg] = Calc_First_order_statics(img_data, mask, bins)
% By Bumwoo Park
% Update: 2018-01-23
% E-mail: julius0628@gmail.com

stat = [];
status = 1;
msg = '';

% 0. initial checks
if sum(size(img_data) ~= size(mask))
    status = -1;
    msg = 'image and mask have different matrix size';
    return;
end

% 1. Extract data in Mask
data_cur = double(img_data(mask));
Counts = size(data_cur, 1);
Max = max(data_cur);
Min = min(data_cur);
Range = Max - Min;

[nelements,centers] = hist(data_cur, bins);
prob = nelements / Counts;
prob = prob(prob>eps); % exclude zero

% 2. 
stat.Range_cover = Range / (max(img_data(:)) - min(img_data(:)) + eps);
stat.Energy = sumsqr(data_cur);
stat.Entropy = -sum(prob .* log(prob + eps) ./ log(2));
stat.Kurtosis = kurtosis(data_cur) - 3;
stat.Max = max(data_cur);
stat.Mean = mean(data_cur);
stat.Mad = mad(data_cur);
stat.Median = median(data_cur);
stat.Min = min(data_cur);
stat.Counts = size(data_cur, 1);
stat.Range = Range;
stat.Rms = rms(data_cur);
stat.Skewness = skewness(data_cur);
stat.Std = std(data_cur);
stat.Sum = sum(data_cur);
stat.Uniformity = sum(prob.*prob);
stat.Var = var(data_cur);
end