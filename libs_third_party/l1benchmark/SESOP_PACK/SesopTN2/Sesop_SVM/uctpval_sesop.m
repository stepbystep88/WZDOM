function  [f,g,hP] = uctpval_sesop(x,P,par)

par.pt=1;  % Reformulate to scalar problem
eps_hess=1e-8;

if nargout==1, [f] = uctpval(x,par);end
if nargout>=2, [f,g] = uctpval(x,par);end

if nargout>=3,
	[np,mp]=size(P);
	hP=zeros(np,mp);
	for ihess=1:mp
		x1=x+eps_hess*P(:,ihess);
		[f2,g2]=uctpval(x1,par);
		hP(:,ihess)=(1/eps_hess)*(g2-g);
	end
end