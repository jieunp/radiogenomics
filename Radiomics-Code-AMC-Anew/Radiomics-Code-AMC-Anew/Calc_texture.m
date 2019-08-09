function [texture, status, msg] = Calc_texture(img_data, mask, bins, max_dist, mask_threshold, mask_res, disp_prog, type, alpha)
% By Bumwoo Park
% Update: 2018-01-23
% E-mail: julius0628@gmail.com


status = 1;
msg = '';

%% For None-wavelet filter
cur_wavelet_txt = 'NONE';
% 1. Calc First-order statics
if disp_prog
    disp('For None-wavelet........');
    disp('Calculating First-order statics....');
end
[stat_FO, status, msg] = Calc_First_order_statics(img_data, mask, bins);

% 2. Calc Shape Features
if disp_prog
    disp('Calculating Shape Features....');
end
[stat_Shape, status, msg] = Calc_Shape(mask, mask_res);

% 3. Calc Run-Length features
if disp_prog
    disp('Calculating Run-Length Features....');
end
[stat_RL, status, msg] = Calc_RL(img_data, mask, bins, type, alpha);

% 4. Calc GLCM features by distance
for dist = 1 : max_dist
    if disp_prog
        disp(['Calculating GLCM Features (distance = ' num2str(dist) ')....']);
    end
    [stat_GLCM{dist}, status, msg] = Calc_GLCM(img_data, mask, bins, dist, type, alpha);
end

texture.(cur_wavelet_txt).FO = stat_FO;
texture.(cur_wavelet_txt).SHAPE = stat_Shape;
texture.(cur_wavelet_txt).GLCM = stat_GLCM;
texture.(cur_wavelet_txt).RL = stat_RL;

%% For wavelet filter
if ndims(img_data) == 3
    % 1. Apply wavelet transform
    wavelet_img = wavedec3(img_data, 1, 'Haar');
    % Ordering
    wavelet_String = {'LLL','HLL','LHL','HHL','LLH','HLH','LHH','HHH'};
    wavelet_Order = [1 5 3 7 2 6 4 8];
    
    % 2. interpolate mask as half size
    [y, x ,z]= ndgrid(linspace(1,size(mask,1),round(size(mask,1)/2)), linspace(1,size(mask,2),round(size(mask,2)/2)), linspace(1,size(mask,3),round(size(mask,3)/2)));mask_half_tmp=interp3(single(mask),x,y,z);
    if mask_threshold ~= 0
        mask_half = mask_half_tmp>=mask_threshold;
    else
        mask_half = mask_half_tmp>mask_threshold;
    end
    
    % 3. Calc Each texture
    for itr = 1 : length(wavelet_img.dec)
        cur_idx = wavelet_Order(itr);
        if disp_prog
            disp(['For Wavelet (' wavelet_String{cur_idx} ')....']);
        end
        cur_img = wavelet_img.dec{cur_idx};
        cur_wavelet_txt = wavelet_String{cur_idx};
        
        % 1. Calc First-order statics
        if disp_prog
            disp('Calculating First-order statics....');
        end
        [stat_FO, status, msg] = Calc_First_order_statics(cur_img, mask_half, bins);
        
        % 2. Calc Run-Length features
        if disp_prog
            disp('Calculating Run-Length Features....');
        end
        [stat_RL, status, msg] = Calc_RL(cur_img, mask_half, bins, type, alpha);
        
        % 3. Calc GLCM features by distance
        for dist = 1 : max_dist
            if disp_prog
                disp(['Calculating GLCM Features (distance = ' num2str(dist) ')....']);
            end
            [stat_GLCM{dist}, status, msg] = Calc_GLCM(cur_img, mask_half, bins, dist, type, alpha);
        end
        
        texture.(cur_wavelet_txt).FO = stat_FO;
        texture.(cur_wavelet_txt).GLCM = stat_GLCM;
        texture.(cur_wavelet_txt).RL = stat_RL;
    end
end
end