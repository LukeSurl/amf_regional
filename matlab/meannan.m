%Code to compute mean for a multidimensional matrix that contains missing data
%Modified from mean.m and nanmean.m as written by MATLAB

%Randall Martin
%April 10, 1999

function y = meannan(x,dim)

if isempty(x) 			% Check for empty input.
    y = NaN;
    return
end

if nargin==1, 
  % Determine which dimension SUM will use
  dim = min(find(size(x)~=1));
  if isempty(dim), dim = 1; end
end

% Replace NaNs with zeros.
nans = isnan(x);
i = find(nans);
x(i) = zeros(size(i));

if min(size(x))==1,
  count = length(x)-sum(nans,dim);
else
  count = size(x,dim)-sum(nans,dim);
end

% Protect against a column of all NaNs
i = find(count==0);
count(i) = ones(size(i));
y = sum(x,dim)./count;
y(i) = i + NaN;
