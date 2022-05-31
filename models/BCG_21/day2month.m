function m = day2month(d)


ncols = size(d,2);

nrows = size(d,1);

remainder = mod(nrows,30);
if remainder ~=0
    addrows = 30-remainder;
    
    d = [d;nan(addrows,ncols)];

    nrows = size(d,1);

end

m = nan(nrows/30,ncols);

for thiscol=1:ncols
    m(:,thiscol) = transpose(nanmean(reshape(d(:,thiscol),30,nrows/30)));
    
end