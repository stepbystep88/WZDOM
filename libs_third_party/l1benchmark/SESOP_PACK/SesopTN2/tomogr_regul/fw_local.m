function Hw_local = fw_local(w_local,par)    %indicator_wind,lam,mu,angles,res)
global Niter
w=zeros(par.Nbins,par.Nang);
for i=1:par.Nang, w(par.iiw(:,i),i)=w_local(:,i);end


Atw=myadjradon(w,par);
A2A2tw=myradon((1-par.wind).*Atw,par);
Hw=A2A2tw + par.lam*w + par.b*(par.b(:)' * w(:));

for i=1:par.Nang, Hw_local(:,i)= Hw(par.iiw(:,i),i);end

