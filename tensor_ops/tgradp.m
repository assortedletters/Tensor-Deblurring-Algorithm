%Performs a forwards finite difference gradient approximation
function [hor vert] = tgradp(x)
[m n p] = size(x);
hor = zeros(m,n,p);
vert = zeros(m,n,p);

hor(:,1:end-1,:) = x(:,2:end,:) - x(:,1:end-1,:);
hor(:,end,:) = x(:,1,:) - x(:,end,:);
vert(1:end-1,:,:) = x(2:end,:,:) - x(1:end-1,:,:);
vert(end,:,:) = x(1,:,:) - x(end,:,:);
