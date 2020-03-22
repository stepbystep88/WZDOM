function [d] = quasi_newt_mz(g,DG,DX,Hdiag)
% BFGS Search Direction
%
% This function returns the (L-BFGS) approximate inverse Hessian,
% multiplied by the gradient
%
% If you pass in all previous directions/sizes, it will be the same as full BFGS
% If you truncate to the k most recent directions/sizes, it will be L-BFGS
%
% s - previous dg search directions (p by k)
% y - previous dx step sizes (p by k)
% g - gradient (p by 1)
% Hdiag - value of initial Hessian diagonal elements (scalar)




g_alpha=(g'*DX)'
h_alpha= DX'*DG;  % Hessian in alpha

%%% Solve modified Newton system: Invert sign of negative eigenvalues %%%%%

[S,Lambda]=eig(h_alpha);
lam=diag(Lambda);
lam=abs(lam);
lam=max(lam,1e-12*max(lam));
d_alpha=-S*(diag(1./lam)*(S'*g_alpha)); %Newton direction in alpha

