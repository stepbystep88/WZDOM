function H=hadamard_matrix(N)
%Create Hadamard matrix NxN

H0=[1 1 ; -1 1]/sqrt(2);
H=H0;
for i=1:1:log(N)/log(2)-1,
  H=kron(H,H0);
end;