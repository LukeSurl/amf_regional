%crude lowpass filter

function data = smooth(noisy,nstep)

ln = length(noisy);
half = floor(nstep/2);
for i = 1:nstep:ln-nstep
  ind = 1+floor(i/nstep);
  val(ind) = i;
  mn(ind) = nanmean(noisy(i:i+nstep-1));
end
data = interp1(half+val,mn,1:ln,'spline');
   
