function [I] = teye(height,depth)
I = zeros(height,height,depth);
I(:,:,1) = eye(height);
end