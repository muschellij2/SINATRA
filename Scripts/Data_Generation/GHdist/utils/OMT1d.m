function OMTDist = OMT1d(mu,nu,p)
%OMT1D 1-dimensional Wasserstein distance
%
% Tingran Gao (tingrangao@galton.uchicago.edu)
% last modified: Oct 28, 2017
%

if (nargin < 3)
    p = 1;
end

if (p == 1)
    %%% In dimension 1, Wasserstein-1 distance is the same as L1 distance
    s = mu; s = s/sum(s);
    t = nu; t = t/sum(t);        
    OMTDist = abs(cumsum(s)-cumsum(t));
else
    lattice = 0:0.01:1;
    %%% In dimension 1, Wasserstein-p (p>1) distance is the L1 distance
    %%% between inverses of cumulative distribution functions
    s = mu; s = s/sum(s);
    t = nu; t = t/sum(t);
    cdf_s = cumsum(s);
    cdf_t = cumsum(t);
    %         cdf_s = interp1(1:length(cumsum(s)),cumsum(s),1:0.5:length(cumsum(s)));
    %         cdf_t = interp1(1:length(cumsum(t)),cumsum(t),1:0.5:length(cumsum(t)));
    is = invCDF(cdf_s, lattice);
    it = invCDF(cdf_t, lattice);
    OMTDist = (sum(abs(is-it).^p))^(1/p);
end

end


function is = invCDF(s,l)
%%% s is a cdf defined equi-distantly on [1,2,\cdots,length(s)]
%%% l is strictly increasing and equi-distant points from 0 to 1
%%% l(1) = 0 and l(end) = 1
if (l(1) ~= 0)
    error('second argument should start with 0');
end
if (l(end) ~= 1)
    error('second argument should end with 1');
end
is = zeros(1,length(l)-1);
for j=1:(length(l)-1)
    is(j) = find(s>l(j), 1)/length(s);
end
end

