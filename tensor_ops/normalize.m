function [V a] = normalize(X)
X = double(X);
tol = 1e-5;
dims = size(X);
if (dims(2)~=1 || length(dims) < 3)
    error('not a twisted matrix');
end
a = zeros(1,1,dims(3));
V = fft(X,[],3);
for j = 1:dims(3)
    a(j) = norm(V(:,:,j),2);
    if a(j) > tol
        V(:,:,j) = V(:,:,j)/a(j);
    else
        %giving it junk
        V(:,:,j) = randn(dims(3),1); a(j) = norm(V(:,:,j),2);
        V(:,:,j) = V(:,:,j)/a(j); a(j) = 0;
    end
end
V = ifft(V,[],3); a = ifft(a,[],3);