function [stat, status, msg] = Calc_RL(img_data, mask, bins, type, alpha)
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
max_RL = max(size(img_data));
for idx_dir = 1 : length(Dir_X)
   Mat_cur =  zeros(bins, max_RL);
   Mat_visit = zeros(size(img_data));   
   
   % 5. Loop for each voxels
   for idx_voxel = 1 : Ele_Mask
       
       % 5.1 Check visited voxel
       if (Mat_visit(X_idx(idx_voxel), Y_idx(idx_voxel), Z_idx(idx_voxel)) == 1)
          continue; 
       end
       
       % 5.2 Get current signal intensity           
       SI_cur = img_quantized(X_idx(idx_voxel), Y_idx(idx_voxel), Z_idx(idx_voxel));
       cur_length = 1;
       for dir_sign = [1 -1]
           for dist = 1 : max_RL
               % 5.3 Check within image boundary
               Next_X = X_idx(idx_voxel) + dir_sign * dist * Dir_X(idx_dir);
               Next_Y = Y_idx(idx_voxel) + dir_sign * dist * Dir_Y(idx_dir);
               Next_Z = Z_idx(idx_voxel) + dir_sign * dist * Dir_Z(idx_dir);
               if Next_X < 1 || Next_X > size(img_data, 1)
                   continue;
               end
               if Next_Y < 1 || Next_Y > size(img_data, 2)
                   continue;
               end
               if Next_Z < 1 || Next_Z > size(img_data, 3)
                   continue;
               end
               
               % 5.4 Check outside mask
               if ~mask(Next_X, Next_Y, Next_Z)
                   break;
               end
               
               % 5.5 Get next signal intensity
               SI_next = img_quantized(Next_X, Next_Y, Next_Z);
               
               % 5.6 Compare current and next signal intensity
               if SI_cur == SI_next
                   cur_length = cur_length+1;
                   Mat_visit(Next_X, Next_Y, Next_Z) = 1;
               else
                   break;
               end
           end
       end
       
       % 5.3 fill GLCM matrix
       Mat_cur(SI_cur, cur_length) =  Mat_cur(SI_cur, cur_length) + 1;       
   end  
   
   % 6. Add matrix along with all directions
   Mat_RL{idx_dir} = Mat_cur';
end

% 7. Calc RL for each directions
stat_tmp = zeros(length(Mat_RL), 12);
for itr = 1 : length(Mat_RL)
    [stat_tmp(itr,:), tGLRLM] = GLRLM_Features1(Mat_RL(itr));
end

% 8. return features..
stat.GLN = mean(stat_tmp(:,3));
stat.GLN_std = std(stat_tmp(:,3));

stat.HGRE = mean(stat_tmp(:,7));
stat.HGRE_std = std(stat_tmp(:,7));

stat.LRE = mean(stat_tmp(:,2));
stat.LRE_std = std(stat_tmp(:,2));

stat.LRHGE = mean(stat_tmp(:,11));
stat.LRHGE_std = std(stat_tmp(:,11));

stat.LRLGE = mean(stat_tmp(:,10));
stat.LRLGE_std = std(stat_tmp(:,10));

stat.LGRE = mean(stat_tmp(:,6));
stat.LGRE_std = std(stat_tmp(:,6));

stat.NRUN = mean(stat_tmp(:,12));
stat.NRUN_std = std(stat_tmp(:,12));

stat.RLN = mean(stat_tmp(:,4));
stat.RLN_std = std(stat_tmp(:,4));

stat.RP = mean(stat_tmp(:,5));
stat.RP_std = std(stat_tmp(:,5));

stat.SRE = mean(stat_tmp(:,1));
stat.SRE_std = std(stat_tmp(:,1));

stat.SRHGE = mean(stat_tmp(:,9));
stat.SRHGE_std = std(stat_tmp(:,9));

stat.SRLGE = mean(stat_tmp(:,8));
stat.SRLGE_std = std(stat_tmp(:,8));

end