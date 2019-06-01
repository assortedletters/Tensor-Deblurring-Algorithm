% performs tensor transpose operation on a 3rd order tensor.
% transposes each frontal slice, and reverses the order of slices from
% second frontal slice to the final frontal slice. 
function X = ttran(x)
dims = size(x);
x = double(x);
X = zeros([dims(2),dims(1),dims(3)]);
X(:,:,1) = x(:,:,1)';
k = 2;
for i = dims(3):-1:2
    X(:,:,k) = x(:,:,i)';
    k = k+1;
end
X = (X);