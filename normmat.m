function res = normmat(x)
% NORMMAT returns the euclidian norm of a matrix.
%
%   x     Input matrix
%   res   Euclidian norm of x

res = sqrt(sum(sum(x.^2)));

end

