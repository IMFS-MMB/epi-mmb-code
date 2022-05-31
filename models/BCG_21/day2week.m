function w = day2week(d)


ncols = size(d,2);

nrows = size(d,1);

remainder = mod(nrows,7);
if remainder ~=0
    addrows = 7-remainder;
    
    d = [d;nan(addrows,ncols)];

    nrows = size(d,1);

end

w = nan(nrows/7,ncols);

for thiscol=1:ncols
    w(:,thiscol) = transpose(nanmean(reshape(d(:,thiscol),7,nrows/7)));
    
end