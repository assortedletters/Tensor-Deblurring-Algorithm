function V = Videorgb2gray(v)
[m n p q] = size(v);
V = zeros(m,n,1,q);
if (p~=1)
    for i = 1:q
        V(:,:,1,i) = rgb2gray(v(:,:,:,i));
    end
else
    V = v;
end