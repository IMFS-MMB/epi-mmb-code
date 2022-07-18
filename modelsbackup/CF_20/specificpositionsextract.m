function out = specificpositionsextract(input)

% this part is written by Kailong Liu to extract the relevent results
tmp1 = length(input)/21;
tmp2 = ([1:tmp1]-1)*21 + 1;
out = input(tmp2);