function q = month2quarter(m)


ncols = size(m,2);

nrows = size(m,1);

remainder = mod(nrows,3);
if remainder ~=0
    addrows = 3-remainder;
    
    m = [m;nan(addrows,ncols)];

    nrows = size(m,1);

end

q = nan(nrows/3,ncols);

for thiscol=1:ncols
    q(:,thiscol) = transpose(nanmean(reshape(m(:,thiscol),3,nrows/3)));
    
end