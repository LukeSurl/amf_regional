function data = smooth_new(noisy,nanlat)

newlat = nanlat(~isnan(noisy));
noisy(isnan(noisy)) = [];

sizenoisy = size(noisy,2);
data = NaN.*ones(size(nanlat));

for i = 1:length(nanlat)
    
    [a lati] = min(abs(newlat-nanlat(i)));
    data(i) = mean(noisy(max([1 lati-10]):min([lati+10 sizenoisy])));
    
end
