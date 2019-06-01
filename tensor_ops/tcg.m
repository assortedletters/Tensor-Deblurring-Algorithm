% Assume tensor A is positive definite blurring kernel,
% matrix B is image to be deblurred
function [X a] = tcg(A,B)
% if not twisted, we twist!
if(size(B,2) ~= 1)
    B = twist(B);
end
X = zeros(size(B));
[R{1} a] = normalize(B); P{1} = R{1};
for i = 1:size(A,3)
    left = tinverse(tprod(tprod(ttran(P{i}),A),P{i}));
    right = tprod(ttran(R{i}),R{i});
    c = tprod(left,right);
    X(:,i+1,:) = X(:,i,:) + tprod(P{1},c);
    R{i+1} = R{i} - tprod(A,tprod(P{1},c));
    left = tinverse(tprod(ttran(R{i}),R{i}));
    right = tprod(ttran(R{i+1}),R{i+1});
    d = tprod(left,right);
    P{i+1} = R{i+1} - tprod(P{i},d);
end
X = tprod(X(:,end,:),a);