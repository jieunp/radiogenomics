function [stat, status, msg] = Calc_GLCM(img_data, mask, bins, dist, type, alpha)
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

% 1. Check dimension of image
dim_num = ndims(img_data);
if dim_num == 2
    Dir_X = [1 1 0 -1];
    Dir_Y = [0 1 1 1];
    Dir_Z = [0 0 0 0];
elseif dim_num == 3
    Dir_X = [0 -1 -1 -1 0 0 0 -1 1 -1 1 -1 1];
    Dir_Y = [1 1 0 -1 1 0 -1 0 0 1 -1 -1 1];
    Dir_Z = [0 0 0 0 -1 1 -1 -1 -1 -1 -1 -1 -1];
else
    status = -1;
    msg = 'image and mask are not 2D or 3D';
    return;
end

% 2. Quantize image
[img_quantized, status, msg] = Quantized_img(img_data, mask, bins, type, alpha);

% 3. Extract image data on Mask.
On_Mask = find(mask~=0);
Ele_Mask = size(On_Mask, 1);
[X_idx, Y_idx, Z_idx] = ind2sub(size(img_data), On_Mask);

% 4. Loop for each directions
Mat_GLCM = zeros(bins, bins, length(Dir_X));
for idx_dir = 1 : length(Dir_X)
   Mat_cur =  zeros(bins, bins);
   
   % 5. Loop for each voxels
   for idx_voxel = 1 : Ele_Mask
       % 5.1 Check within image boundary
       Next_X = X_idx(idx_voxel) + dist * Dir_X(idx_dir);
       Next_Y = Y_idx(idx_voxel) + dist * Dir_Y(idx_dir);
       Next_Z = Z_idx(idx_voxel) + dist * Dir_Z(idx_dir);
       if Next_X < 1 || Next_X > size(img_data, 1)
           continue;
       end
       if Next_Y < 1 || Next_Y > size(img_data, 2)
           continue;
       end
       if Next_Z < 1 || Next_Z > size(img_data, 3)
           continue;
       end
       
       % 5.2 Get current and next signal intensity
       SI_cur = img_quantized(X_idx(idx_voxel), Y_idx(idx_voxel), Z_idx(idx_voxel)); 
       SI_next = img_quantized(Next_X, Next_Y, Next_Z);
       
       % 5.3 fill GLCM matrix
       Mat_cur(SI_cur, SI_next) =  Mat_cur(SI_cur, SI_next) + 1;
       Mat_cur(SI_next, SI_cur) =  Mat_cur(SI_next, SI_cur) + 1;       
   end  
   
   % 6. Add matrix along with all directions
   Mat_GLCM(:, :, idx_dir) = Mat_cur;
end

% 7. Calc GLCM for each direction
for itr = 1 : length(Dir_X)
    stat_total(itr) = GLCM_Features3(Mat_GLCM(:, :, itr), 0);
end

% 8. Fill mean and std of each features...
names = fieldnames(stat_total);
for itr = 1 : length(names)
    cur_arr = [stat_total.([names{itr}])];
    stat.([names{itr}]) = mean(cur_arr);
    stat.([names{itr} '_std']) = std(cur_arr);    
end

end