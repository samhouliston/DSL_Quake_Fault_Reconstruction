function res = gaussianValue(X, m, covar)
% Computes probability of dataset X (each _column_ is one datum),
% given a gaussian pdf, Normal(m, covar).
%
%   Example:
%           gaussianValue([4 3 2; 2 3 4], [3 3]', eye(2))
%       will return [0.0585 0.1592 0.0585]'.
%
% G.Sfikas 7 feb 2007.
[d,N]   = size(X); 
detC    = det(covar);
if(detC<=0)
    res     = zeros(N,1);
    return;
end
    bigM = (X-m * ones(1, N))';
    %mah = sum((X-bigM)' * (inv(covar)+eps*eye(size(X,1))) .* (X-bigM)', 2);
    %mah = sum(bigM * inv(covar) .* bigM, 2);
    mah = sum(bigM * (inv(covar)+eps*eye(d)) .* bigM, 2);
    res = (2*pi)^(-0.5*d) * detC^(-0.5) * exp(-0.5*mah); 

end