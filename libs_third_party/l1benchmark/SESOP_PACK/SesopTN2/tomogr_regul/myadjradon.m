% Adjoint Radon Transform

function  x = myadjradon(y,par)
  
if par.flag_fastradon
	x = fastadjradon(y,par);
else
    y=y';
    xtest=fastbp(y(1,:),1,4);
	[n,n]=size(xtest);
	x=zeros(n);

	 for i = 1:length(par.angles)
        angle = par.angles(i);
		x = x + fastbp(y(i,:),angle,4);
	 end
end

	
	
