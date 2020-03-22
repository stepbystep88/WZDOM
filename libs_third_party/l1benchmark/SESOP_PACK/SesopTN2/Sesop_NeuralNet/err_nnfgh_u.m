function [r,g_r,HrZ]=err_nnfgh_u(vcub,Z,par,y)
% Error function of Single-layer feed-forward Neural Net as a function of to v, c, U=WX, b;
% its gradient and Hessian-vector (matrix) multiplication
%
% Call: [r,g_r,HrZ]=err_nnfgh_u(vcub,Z,par,y)
%
% Input:
%   vcub - v,c,U,b in a single vector, where U=WX; 
%   Z - matrix to be multiplied by Hessian (if needed)
%   par - user parameters
%   y - vector of desired NN responce to X
%   
% Output:
%    r   -  NN Error: 1/2||ynn -y||^2, where ynn=v'Phi(U+b*1')+c;  
%    g_r -  gradient of r with respevt to vcub
%    HrZ -  Hessian*Z  

% Michael Zibulevsky, 24.07.2008; 29.07.2008; 05.08.2008
%
% Copyright (c) 2008. All rights reserved. Free for academic use. No warranty 


M=par.nneurons;
K=par.Ktrain_samples;

indv=[1:M];
indc=M+1;
%indU=[M+2: M+1+M*K];
indb=[(M*(K+1) +2) : (M*(K+2) +1)];

v=vcub(indv);
c=vcub(indc);
u=vcub(M+2: M+1+M*K);
U=reshape(u, M, K);
b=vcub(indb);

ones1K=ones(1,K);
onesK1=ones(K,1);

U=U+b*ones1K;  % Add bias


if nargout==1
	[Phi]=sigmoid_mz(U, par.eps_sigmoid, par.flagXnew);
elseif nargout==2
	[Phi,Phi1]=sigmoid_mz(U, par.eps_sigmoid, par.flagXnew);
elseif nargout==3
	[Phi,Phi1,Phi2]=sigmoid_mz(U, par.eps_sigmoid, par.flagXnew);
end


F= (v'*Phi)' +c - y(:);
r=0.5*sumsqr(F);

if nargout >1,     %Gradients
	vFt=v*F';
	gr_v=Phi*F;  
	gr_c=sum(F); 
	Gr_U=vFt.*Phi1;
	gr_b=Gr_U*onesK1;
	g_r=[gr_v;gr_c; Gr_U(:); gr_b];
end

if nargout >2, % Multiply Hessian by a vector (matrix) Z
	ones2Mp1=ones(2*M+1,1);
	J=size(Z,2);
	HrZ=zeros(size(Z));
	for j=1:J
		H1z=0;
		z_v=Z(indv,j);
		z_c=Z(indc,j);
		z_b=Z(indb,j);
		z_u=Z(M+2: M+1+M*K,j);
		Z_u=reshape(z_u,size(U)) +z_b*ones1K;
		
		Phi1Z_u=Phi1.*Z_u;
		dF= (   z_v'*Phi + ones1K*z_c + v'*Phi1Z_u  )';
		Hz_v=Phi1Z_u*F + Phi*dF;
		Hz_c=sum(dF);
		HZ_u= (z_v*F' + v*dF').*Phi1 + vFt.*Phi2 .*Z_u;
		HrZ(1:M,j)=Hz_v;
		HrZ(M+1,j)=Hz_c;
		HrZ(M+2:end-M,j)=HZ_u(:);
		HrZ(indb,j)=HZ_u*onesK1;
	end
end




