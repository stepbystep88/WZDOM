function report(w, iter,par)

%if mod(iter,5)==0,
	Atw=myadjradon(w,par);
	%figure(1)
	subplot(121);imagesc(w); colorbar;title(sprintf('iteration %d', iter));
	subplot(122);imagesc((Atw)); colorbar;
	title(sprintf('NormW = %g  normA2tw=%g lam=%g',norm(w(:))/(par.wind(:)'*Atw(:)),norm(Atw(:).*(1-par.wind(:))) /(par.wind(:)'*Atw(:)),par.lam));
    drawnow
%end
