%converts an m x 1 x n tensor into an m x n matrix. The inverse of twist
function x = tsqueeze(X)
dims = size(X);
if (dims(2) ~= 1 || length(dims) < 3)
    error('not a twisted matrix');
end
x = zeros(dims(1),dims(3));
for i = 1:dims(3)
    x(:,i) = X(:,1,i);
end
