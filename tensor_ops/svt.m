function W = svt(X,tau)
X = double(X);
if(length(size(X)) == 3)
    W = fft(X,[],3);
    n = size(W,3);
    for i = 1:n
        [u s v] = svd(W(:,:,i));
        s = s-tau*eye(size(s));
        s(s<0) = 0;
        U(:,:,i) = u; S(:,:,i) = s; V(:,:,i) = v;
        
    end
    U = ifft(U,[],3); S = ifft(S,[],3); V = ifft(V,[],3);

end
if(length(size(X)) == 4)
    n = size(W,4);
    W = fft(X,[],4);

    for i = 1:n
        [u s v] = svd(W(:,:,:,i));
        s = s-tau*eye(size(s));
        s(s<0) = 0;
        U(:,:,:,i) = u; S(:,:,:,i) = s; V(:,:,i) = v;
    end
    U = ifft(U,[],4); S = ifft(S,[],4); V = ifft(V,[],4);

end
W = tprod(tprod(U,S),ttran(V));
W = double(W);
