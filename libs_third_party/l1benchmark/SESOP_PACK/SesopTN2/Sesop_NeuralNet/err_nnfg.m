function [f,g]=err_nnfg(vcwb,par,X,y)
%Error function of Single-layer feed-forward Neural Net and its gradient with respect to weights   
%
%Call:  [f,g]=err_nnfg(vcwb,par)
%
%Let ynn=v'Phi(WX+b*1'); err_nn=1/2||ynn -y||^2;


% Michael Zibulevsky  24.07.2008 
%
% Copyright (c) 2008. All rights reserved. Free for academic use. No warranty 

par.flagXnew=1;

vcub=multWX(vcwb,par,X);

if nargout==1,
	f=err_nnfgh_u(vcub, [], par, y);
else
	[f, g_u]=err_nnfgh_u(vcub, [], par, y);
	g=multUXt(g_u,par,X);
end




%%%%%%% Old version %%%%%%%%%%%

% [N, K]=size(X);            % matrix X:  N inputs (including ones) x K examples
% 
% M=par.nneurons;
% v=vcwb(1:M);
% c=vcwb(M+1);
% w=vcwb(M+2 : end-M);
% b=vcwb(end-M+1:end);
% 
% W=reshape(w,par.nneurons,N);
% 
% 
% [Phi,Dphi]=sigmoid_mz(W*X+b*ones(1,K), par.eps_sigmoid);
% 
% r=v' * Phi +c - y;  %residual
% 
% f=0.5*sumsqr(r) + 0.5*par.quadrpenpar*sumsqr(vcw);  
% 
% g_v=Phi*r';
% 
% g_c=sum(r(:));
% 
% g_W=zeros(size(W));
% g_b=zeros(M,1);
% 
% %Gr_U=vFt.*Phi1;
% %gr_b=Gr_U*onesK1;
% 
% for i=1:K,
% 	tmp=r(i)*(Dphi(:,i).*v);
%    g_W=g_W+ tmp*X(:,i)';
% end
% 
% g=[g_v; g_c; g_W(:)] + par.quadrpenpar*vcw;
% 

