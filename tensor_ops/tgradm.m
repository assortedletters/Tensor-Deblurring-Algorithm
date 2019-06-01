%Performs a backwards finite difference gradient approximation

function [hor vert] = tgradm(x)
[m n p] = size(x);
hor = zeros(m,n,p);
vert = zeros(m,n,p);

hor(:,2:end,:) = x(:,2:end,:) - x(:,1:end-1,:);
hor(:,1,:) = x(:,1,:) - x(:,end,:);
vert(2:end,:,:) = x(2:end,:,:) - x(1:end-1,:,:);
vert(1,:,:) = x(1,:,:) - x(end,:,:);
