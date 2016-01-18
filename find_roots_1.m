% In the calculation of the spectrum, we need the roots of the zeroth
% Bessel function of the first kind. This program is designed to calculate
% the roots of this function.  Here we calculate the roots of the Bessel
% function of the first kind.


maxs = 500;
roots_1 = zeros(1, maxs);
delta = zeros(1, maxs-1);
incr = pi;

roots_1(1,1) = fzero(@(x)besselj(1,x),2.8);

for s=2:maxs
    roots_1(1, s) = fzero(@(x)besselj(1,x),[roots_1(1,s-1)+incr-0.3, roots_1(1,s-1)+incr+0.3]);
    delta(1, s-1) = roots_1(1, s)-roots_1(1, s-1);
end    

save roots_1.mat roots_1;