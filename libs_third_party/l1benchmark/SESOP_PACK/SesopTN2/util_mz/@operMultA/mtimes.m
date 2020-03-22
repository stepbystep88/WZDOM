function res = mtimes(A,x)

if A.adjoint == 0 %A*x
    res = A.par.multA(x,A.par);
else %At*x
   res = A.par.multAt(x,A.par);
end
