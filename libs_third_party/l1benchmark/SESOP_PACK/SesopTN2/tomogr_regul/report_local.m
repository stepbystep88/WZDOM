function report_local(w_local, iter,par)

w=zeros(par.Nbins,par.Nang);
for i=1:par.Nang, w(par.iiw(:,i),i)=w_local(:,i);end

%if mod(iter,5)==0,
Atw=myadjradon(w,par);

tmp=par.wind(:)'*Atw(:); % normalization 
Atw=Atw/tmp; w=w/tmp; 

%figure(1)
subplot(121);imagesc(w); colorbar;title(sprintf('W: iter= %d', iter));
subplot(122);imagesc((Atw)); colorbar;
title(sprintf('NormW=%g  normA2tw=%g lam=%g', norm(w(:)),norm(Atw(:).*(1-par.wind(:))),par.lam));
drawnow
%end
