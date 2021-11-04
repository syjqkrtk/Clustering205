%%%%% Distance matrix display
%%%%% Byung Chang Chung, KAIST UMLS
%%%%%
%%%%% input: symmetric distance matrix (dmat)
%%%%% output: plotted result
%%%%%

function [a] = visual_dmat(dmat, cn)

[row col] = size(dmat);

max_dist = max(max(dmat));
norm_dmat = uint8(dmat./max_dist.*255);
norm_dmat3 = uint8(zeros(row,col,3));
norm_dmat3(:,:,1) = norm_dmat;
norm_dmat3(:,:,2) = norm_dmat;
norm_dmat3(:,:,3) = norm_dmat;

image(norm_dmat3);
title(cn);
axis image;

end