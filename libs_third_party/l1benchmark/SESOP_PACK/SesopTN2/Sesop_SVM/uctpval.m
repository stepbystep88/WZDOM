function  [F, dF] = uctpval(x,par)
% Unconstrained minimization test problems
% Hans Bruun Nielsen,  IMM, DTU.  00.11.24

% Ensure that  x  is a column vector
x = x(:);

if  par.p <= 21    % Least squares problem
  switch  par.p
   case 1    % Linear function.  Full rank
     f = par.A*x - 1;
     if  nargout > 1,  J = par.A; end
   case 2    % Linear function.  Rank 1
     f = par.A*x - 1;
     if  nargout > 1,  J = par.A; end
   case 3    % Linear function.  Rank 1.  Zero cols. and rows
     f = par.A*x - 1;
     if  nargout > 1,  J = par.A; end
   case 4    % Rosenbrock
     f = [10*(x(2) - x(1)^2); 1-x(1)];
     if  nargout > 1,  J = [-20*x(1)  10; -1  0]; end
   case 5    % Helical Valley
     t = atan(x(2)/x(1))/(2*pi);
     if  x(1) < 0,  t = t + .5; end
     nx = norm(x(1:2));
     f = [10*(x(3) - 10*t); 10*(nx - 1); x(3)];
     if  nargout > 1
       J = zeros(3,3);
       K = -50/pi/nx^2;   J(1,:) = [-K*x(2)  K*x(1)  10];
       K = 10/nx;   J(2,1:2) = K*x(1:2)';   J(3,3) = 1;
     end
   case 6    % Powell singular
     s5 = sqrt(5);   s10 = sqrt(10);
     d3 = x(2) - 2*x(3);   d4 = x(1) - x(4);
     f = [x(1)+10*x(2); s5*(x(3) - x(4)); d3^2; s10*d4^2];
     if  nargout > 1
       J = [1 10 0 0; 0 0 s5 -s5
            0 2*d3*[1 -2] 0; 2*s10*d4*[1 0 0 -1]];
     end
   case 7    % Freudenstein and Roth
     x1 = x(1);   x2 = x(2);
     f = [(x1 - x2*(2 - x2*(5 - x2)) - 13)
          (x1 - x2*(14 - x2*(1 + x2)) - 29)];
     if  nargout > 1
       J = [1  (-2 + x2*(10 - 3*x2))
            1  (-14 + x2*(2 + 3*x2))];
     end
   case 8    % Bard
     D = par.uvwy(:,2:3)*x(2:3);
     f = par.uvwy(:,4) - (x(1) + par.uvwy(:,1)./D);
     if  nargout > 1
       F = -par.uvwy(:,1)./D.^2;
       J = -[ones(size(D))  (F*[1 1]).*par.uvwy(:,2:3)];
     end
   case 9    % Kowalik and Osborne
     u = par.yu(:,2);   x1 = x(1);
     N = u.*(u + x(2));   D = u.*(u + x(3)) + x(4);
     f = par.yu(:,1) - x1*N./D;
     if  nargout > 1
       F = -x1*N./D.^2;
       J = -[N./D  x1*u./D  F.*u  F];
     end
   case 10    % Meyer
     D = par.ty(:,1) + x(3);   e = exp(x(2)./D);
     f = x(1)*e - par.ty(:,2);
     if  nargout > 1
       J = [e  x(1)*e./D  -x(1)*x(2)*e./D.^2];
     end
   case 11    % Watson
     g = par.B*x;   x1 = x(1);   x2 = x(2);  n = length(x);
     f = [(par.A*x(2:end) - g.^2 - 1); x1; x2-x1^2-1];
     if  nargout > 1
       J = zeros(31,n);   g = -2*g;
       J(:,1) = [g; 1; -2*x1];   J(31,2) = 1;
       for  j = 2 : n
         J(1:29,j) = par.A(:,j-1) + g.*par.B(:,j);
       end
     end
   case 12    % Box 3-d
     t = par.t;   E = exp(-t*[x(1:2)' 1  10]);
     c = [1; -1; -x(3)*[1; -1]];
     f = E*c;
     if  nargout > 1
       J = [-t.*E(:,1) t.*E(:,2) E(:,3:4)*[-1;1]];
     end
   case 13    % Jennrich and Sampson
     t = par.t;   E = exp(t*x');
     f = 2*(1 + t) - E*[1; 1];
     if  nargout > 1
       J = -(t*ones(1,2)).*E;
     end
   case 14    % Brown and Dennis
     t = par.t;   st = sin(t);
     d1 = x(1) + x(2)*t - exp(t);   d2 = x(3) + x(4)*st - cos(t);
     f = d1.^2 + d2.^2;
     if  nargout > 1
       J = 2*[d1  d1.*t  d2  d2.*st];
     end
   case 15    % Chebyquad
     z = real(acos(2*x' - 1));   f = -par.y;
     m = length(f);   n = length(x);
     for  r = 1 : m,  f(r) = sum(cos(r*z))/n + f(r); end
     if  nargout > 1
       J = zeros(m,n);   d = sqrt(1 - (2*x' - 1).^2);
       for  r = 1 : m
         J(r,:) = 2/n*r*sin(r*z) ./ d;
       end 
     end
   case 16    % Brown almost linear
     n = length(x);   p = prod(x);
     f = [par.A*x-(n+1); p-1];
     if  nargout > 1
       pp = zeros(1,n);
       for  j = 1 : n,  pp(j) = prod(x([1:j-1 j+1:n])); end 
       J = [par.A; pp];  end
   case 17    % Osborne 1
     t = par.ty(:,1);   E = exp(-t*x(4:5)');
     f = par.ty(:,2) - x(1) - E*x(2:3);
     if  nargout > 1
       J = -ones(length(t),5);
       J(:,2:5) = -[E  -x(2)*t.*E(:,1)    -x(3)*t.*E(:,2)];
     end
   case 18    % Exponential fit.  n = 4
     t = par.ty(:,1);   E = exp(t*x(1:2)');
     f = par.ty(:,2) - E*x(3:4);
     if  nargout > 1
       J = -[E  E];
       J(:,1:2) = -[x(3)*t.*E(:,1)   x(4)*t.*E(:,2)];
     end
   case 19    % Exponential fit.  n = 2
     t = par.ty(:,1);   E = exp(t*x');   c = E\par.ty(:,2);
     f = par.ty(:,2) - E*c;
     if  nargout > 1
       A = E'*E;   H = [t t] .* E;
      G = A\ (diag(H'*f) - H'*E*diag(c));
      J = -(E*G + H*diag(c));
     end
   case 20    % Modified Meyer
     D = par.ty(:,1) + x(3);   e = exp(10*x(2)./D - 13);
     f = x(1)*e - par.ty(:,2);
     if  nargout > 1
       q = (10*x(1)./D) .* e;
       J = [e  q  -(x(2)./D).*q];
     end
   case 21    % Separated Meyer
     x2 = par.ty(:,1) + x(2);   a = x(1) ./ x2;
     F = exp(a);   E = F'*F;   c = (F' * par.ty(:,2))/E;
     f = F*c - par.ty(:,2);
     if  nargout > 1
       d1 = x(1) ./ x2;   d2 = c ./ x2;
       by = (f + F*c) ./ x2;   by = [by -by.*d1];   
       dc = -(F'*by)/E;
       J = [F F] .* [dc(1)+d2  dc(2)-d1.*d2];
     end
   end
 
 if  par.pt  % Reformulate to scalar problem
   F = .5 * norm(f)^2;
   if  nargout > 1  % Gradient 
     dF = J'*f;
   end
 else
   F = f;
   if  nargout > 1, dF = J; end  % Jacobian
 end
 
else    % 'Born' scalar function
  switch  par.p
    case 22    % Exp and squares
      e = exp(-sum(x));  ii = (1:length(x)).^2;
      F = e + .5*sum(ii(:) .* x(:).^2);
      if  nargout > 1, dF = -e + ii(:) .* x(:); end
  end
end
