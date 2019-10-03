function [ ob ] = obj( x )
%OBJ Test objective function

% ob = sum((x).^2);
ob = bsAckley(x, 0);




end

