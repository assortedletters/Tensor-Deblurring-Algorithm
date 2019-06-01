function [U S V] = tsvd(A)
dims = size(A);
if (length(dims) ~= 3)
    error('not a 3rd order tensor')
end
A = double(A);
m = dims(1); n = dims(2); p = dims(3);
U = zeros(m,m,p); S = zeros(m,n,p); V = zeros(n,n,p);
A = fft(A,[],3);
for i = 1:dims(3)
    [u s v] = svd(A(:,:,i));
    U(:,:,i) = u; S(:,:,i) = s; V(:,:,i) = v;
end
U = ifft(U,[],3); S = ifft(S,[],3); V = ifft(V,[],3);
U = tensor(U); S = tensor(S); V = tensor(V);