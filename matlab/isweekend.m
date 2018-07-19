ndays = [0,31,30,31,31,30,31,30];    %may, jun, jul, aug, sep, ond
ndays = cumsum(ndays);
ndays2005 = [0,31,28,31,30];
ndays2005 = cumsum(ndays2005);
if eval(idate(5:6)) < 5
  dmonth = ndays2005(eval(idate(5:6)));
else
  dmonth = ndays(eval(idate(5:6))-4);
end
elapsed = eval(idate(7:8)) + dmonth;
if (rem(elapsed,7)==1) | (rem(elapsed,7)==2)
  weekend = true;
else
  weekend = false;
end


