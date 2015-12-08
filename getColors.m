function [dist_colors,flip_color] = getColors
% colormap
N = 64;
CBcols = [0.85 0.85 0.85];%[0.9 0.9 0.9];
t=(1:N) /N;t=t';
cols = (1-t)*CBcols + t*[0.0 0.0 1.0];
cols(cols>1) =1;
dist_colors = cols;
flip_color = [162,20,47]/255';