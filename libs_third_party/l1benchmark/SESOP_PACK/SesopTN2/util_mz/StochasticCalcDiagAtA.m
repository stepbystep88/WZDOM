function  diag_H=StochasticCalcDiagAtA(multAt,size_x,NumTests,params)
% Compute diag A'A: norms of rows of A'  as empirical variance
% of A'n, where n - standard gaussian iid vector
%
%  Zibulevsky

%NumTests=20;
diag_H=0;
for i=1:NumTests, 
  diag_H=diag_H + (multAt(randn(size_x),params)).^2; 
end
diag_H= (1/NumTests) * diag_H;



