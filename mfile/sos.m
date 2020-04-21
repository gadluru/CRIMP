function Im_out = sos(Im_in,dim)

if nargin == 1
    dim = length(size(Im_in));
end
Im_out = sqrt(sum(conj(Im_in).*Im_in,dim));
