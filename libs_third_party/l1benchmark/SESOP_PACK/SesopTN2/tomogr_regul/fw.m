function Hw = fw(w,par)    
Atw=myadjradon(w,par);
A_Dc2_At_w=myradon((1-par.wind).^2 .* Atw,par);
Hw=A_Dc2_At_w + par.lam*w + par.mu*par.b*(par.b(:)' * w(:));

