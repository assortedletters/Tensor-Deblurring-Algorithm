% finds the inverse of an m x m x n tensor
% currently only works for m = 1 i.e. a tubal scalar

function B = tinverse(A)
if (size(A,1) ~= size(A,2))
    error('invalid tensor input. Frontal slices not square');
end
n = size(A,3); m = size(A,2);
A = fft(double(A),[],3);
B = zeros(size(A));
for i = 1:n
    B(:,:,i) = 1./(A(:,:,i));
end
B = ifft(B,[],3);
B = ttran(B);