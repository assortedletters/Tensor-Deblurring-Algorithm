%performs the tensor t-product on tensors A and B, utilizing fourier domain
%for better computational performance. 
function C = tprod(A,B)
dims1 = size(A); dims2 = size(B);
C = 0;
if(length(dims1) == 3)
    if (dims1(2) ~= dims2(1) || dims1(3) ~= dims2(3))
        error('tensor sizes do not agree');
    end
    A = fft(double(A),[],3); B = fft(double(B),[],3);
    C = zeros(dims1(1),dims2(2),dims2(3));
    for i = 1:dims1(3)
        C(:,:,i) = A(:,:,i)*B(:,:,i);
    end
    C = ifft(C,[],3);
    C = tensor(C);
end
if (length(dims1) == 4)
    A = permute(A,[1 2 4 3]); B = permute(B, [1 2 4 3]);
    if (dims1(2) ~= dims2(1) || dims1(3) ~= dims2(3))
        error('tensor sizes do not agree');
    end
    A = fft(double(A),[],3); B = fft(double(B),[],3);
    C = zeros(dims1(1),dims2(2),dims2(3));
    for i = 1:dims1(3)
        C(:,:,i) = A(:,:,i)*B(:,:,i);
    end
    C = ifft(C,[],3);
    C = permute(C,[1 2 4 3]);
    C = tensor(C);
end