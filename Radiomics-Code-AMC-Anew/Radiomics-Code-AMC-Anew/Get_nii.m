function [img_data, img_size, img_res, status, msg] = Get_nii(file_name)
% By Bumwoo Park
% Update: 2018-05-11
% E-mail: julius0628@gmail.com


img_data = [];
img_size = 1;
img_res = '';
status = 1;
msg = '';

% 0. initial checks
if exist(file_name, 'file') ~=  2
    status = -1;
    msg = 'file does not exist';
    return;
end
try
    img_data_struct = load_nii(file_name);
catch err
    img_data_struct = load_untouch_nii(file_name);
end
img_num_dim = img_data_struct.hdr.dime.dim(1);
img_size = img_data_struct.hdr.dime.dim(2:(img_num_dim+1));
img_res = img_data_struct.hdr.dime.pixdim(2:(img_num_dim+1));
img_data = img_data_struct.img;
end