function report_recon(x, iter,par)

%if mod(iter,5)==0,
fprintf('%d ', iter)
subplot(121); imagesc(par.x00);colorbar;title('Phantom');colormap(gray)
subplot(122); imagesc(x);colorbar;title(sprintf('iteration %d StdError=%g', iter, std(x(:)-par.x00(:))));colormap(gray)
	 drawnow
%end
