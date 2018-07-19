% Get the stratospheric column over the Pacific
%
% Assumptions:
%   Zonal invariance of the stratosphere

function [strat] = get_strat_1d(lat,lon,st)

  %go through different regions
  latind = find(lat<-25);
  lonind = find(lon<=280 & lon>=180);
  st1(latind) = meannan(st(latind,lonind),2);
  latind = find(lat<-15 & lat>=-25);
  lonind = find(lon<=250 & lon>=180);
  st1(latind) = meannan(st(latind,lonind),2);
  latind = find(lat<0 & lat>=-15);
  lonind = find(lon<=250 & lon>=160);
  st1(latind) = meannan(st(latind,lonind),2);
  latind = find(lat<15 & lat>=0);
  lonind = find(lon<=250 & lon>=135);
  st1(latind) = meannan(st(latind,lonind),2);
  latind = find(lat<20 & lat>=15);
  lonind = find(lon<=215 & lon>=135);
  st1(latind) = meannan(st(latind,lonind),2);
  latind = find(lat<23 & lat>=20);   %remove Hawaii
  lonind = find((lon<=198 & lon>=135) | (lon<=215 & lon>=203));
  st1(latind) = meannan(st(latind,lonind),2);
  latind = find(lat<38 & lat>=23);
  lonind = find(lon<=215 & lon>=163);
  st1(latind) = meannan(st(latind,lonind),2);
  latind = find(lat<56 & lat>=38);
  lonind = find(lon<=225 & lon>=173);
  st1(latind) = meannan(st(latind,lonind),2);
  latind = find(lat>=56);
  lonind = find(lon<=200 & lon>=173);
  st1(latind) = meannan(st(latind,lonind),2);
  
  strat = NaN.*st1;
  strat(~isnan(st1)) = smooth(st1(~isnan(st1)),5);    %average over 5 degrees in lat
  strat(strat<0)=NaN;


