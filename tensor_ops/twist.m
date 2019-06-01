% converts an m x n matrix into a m x 1 x n tensor
function X = twist(x)
dims = (size(x));
n = length(dims);
if (n > 2)
    if (dims(3) ~= 1)
        error('not a matrix');
    end
end

X = zeros(dims(1),1,dims(2));
if (n < 3)
    for i = 1:dims(1)
        X(i,1,:)= x(i,:);
    end
else
    for i = 1:dims(1)
        X(i,:,:)= x(i,:,:);
    end
end