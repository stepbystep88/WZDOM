function  res = operMultA(par)

res.adjoint = 0;
res.par=par;

% Register this variable as a operMultA class
res = class(res,'operMultA');
