function B = fit_basis_spectra(I,M,method)
% FIT_BASIS_SPECTRA determines the coefficients of the basis spectra in 
% TR-SAXS spectroscopic data modeling. These coefficients obey I~B*M 
% where the approximation is defined in the least squares sense. All the
% elements of B are constrained to be non negative if the chosen method is 
% 'global'.
%
%   I         Spatio-temporal data or transposed temporal singular vectors
%   M         Model vectors
%   method    Method of analysis: 'svd' or 'global'
%
%   B         Matrix of the basis spectra or matrix of coefficients

if strcmp(method,'svd')
    B = mrdivide(I,M);
else
    Nq = size(I,1);
    Rank = size(M,1);
    B = zeros(Nq,Rank);
    for k = 1:Nq
        B(k,:) = lsqnonneg(M',I(k,:)')';
    end
end

end

