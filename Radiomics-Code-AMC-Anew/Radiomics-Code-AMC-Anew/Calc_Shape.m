function [stat, status, msg] = Calc_Shape(mask, mask_res)
% By Bumwoo Park
% Update: 2019-01-23
% E-mail: julius0628@gmail.com

stat = [];
status = 1;
msg = '';

dim_num = ndims(mask);
if dim_num == 2
    stat.Compt1 = -99999;
    stat.Compt2 = -99999;
    stat.Dispro = -99999;
    stat.Sphe = -99999;
    stat.Surface = -99999;
    stat.SVratio = -99999;
    stat.Volume = -99999;
    return;
end

[Surface, labels] = imSurface(mask, mask_res);
not_a_zero_pt = mask>0;
counts = sum(not_a_zero_pt(:));
Volume = counts * prod(mask_res);

stat.Compt1 = Volume / (sqrt(pi) * (Surface ^ (2/3)));
stat.Compt2 = 36 * pi * (Surface^2) / (Volume^3);
stat.Dispro = Surface / ((6 * sqrt(pi) * Volume)^(2/3));
stat.Sphe = (6 * sqrt(pi) * Volume) ^ (2/3) / Surface;
stat.Surface = Surface;
stat.SVratio = Surface / Volume;
stat.Volume = counts * prod(mask_res);
end