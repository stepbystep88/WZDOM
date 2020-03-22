% Discrete Radon transform 

function y = myradon(x,par)

if par.flag_fastradon
	y = fastradon(x,par);
else
	ytest=fastfp(x,1,4);

	y=zeros(length(par.angles),length(ytest));

	for i = 1:length(par.angles)
		angle = par.angles(i);
		y(i,:)= fastfp(x,angle,4);
	end
	y=y';

end


