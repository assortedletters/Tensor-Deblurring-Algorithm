function D = naive_tprod(A,B)
dims1 = size(A); dims2 = size(B);
if (dims1(2) ~= dims2(1) || dims1(3) ~= dims2(3))
    error('tensor sizes do not agree');
end
m = dims1(1); n = dims2(2); p = dims2(3);
C = circulant(A)*tenmat(B,2)';
D = zeros(m,n,p);
for i = 1:p
    D(:,:,i) = C(1+m*(i-1):m*i,1:n);
end
C = tensor(double(C), [m n p]);
D = tensor(D);