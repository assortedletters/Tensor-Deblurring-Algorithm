function S = TNN(X)
X = double(X);
X = fft(X,[],3);
n3 = size(X,3);
S=0;
for i = 1:n3
    [~, s, ~] = svd(X(:,:,i));
    S = S+sum(sum(s));
end
