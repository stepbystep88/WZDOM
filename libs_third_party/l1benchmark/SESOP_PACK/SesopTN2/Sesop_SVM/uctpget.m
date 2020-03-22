function [par, x0, tau0, delta0] = uctpget(pno,m,n)
% Get unconstrained minimization test problem.
% Call:   [par, x0, tau0, delta0] = uctpget(pno,m,n)
% See H.B. Nielsen:  "UCTP - Test Problems for Unconstrained
% Optimization", Report IMM-REP-2000-18

% Hans Bruun Nielsen,  IMM, DTU.  00.11.24

  switch  pno
   case 1    % Linear function.  Full rank
     par = struct('p',1, 'pt',0, 'A',eye(m,n) - (2/m)*ones(m,n), ...
                  'xm',-ones(n,1));
     x0 = ones(1,n);   tau0 = 1e-8;   delta0 = 10;
   case 2    % Linear function.  Rank 1
     xm = (6/(2*m+1)/n/(n+1)) * ones(n,1);
     par = struct('p',2, 'pt',0, 'A',[1:m]'*[1:n], 'xm',xm);
     x0 = ones(1,n);   tau0 = 1e-8;   delta0 = 10;
   case 3    % Linear function.  Rank 1.  Zero cols. and rows
     xm = (6/(2*m-3)/(n-2)/(n+1)) * ones(n,1);
     par = struct('p',3, 'pt',0, 'A',[0 1:m-2 0]'*[0 2:n-1 0], 'xm',xm);
     x0 = ones(1,n);   tau0 = 1e-8;   delta0 = 10;
   case 4    % Rosenbrock
     par = struct('p',4, 'pt',0, 'xm',[1;1]); 
     x0 = [-1.2  1];  tau0 = 1;   delta0 = 1;
   case 5    % Helical Valley
     par = struct('p',5, 'pt',0, 'xm',[1;0;0]); 
     x0 = [-1  0  0];  tau0 = 1;   delta0 = 1;
   case 6    % Powell singular
     par = struct('p',6, 'pt',0, 'xm',zeros(4,1));
     x0 = [3  -1  0  1];   tau0 = 1e-8;   delta0 = 1;
   case 7    % Freudenstein and Roth
     xm = ([53;2] - sqrt(22)*[4;1])/3;
     par = struct('p',7, 'pt',0, 'xm',xm); 
     x0 = [.5  -2];  tau0 = 1;     delta0 = 1;
   case 8    % Bard
     Bard = [0.14  0.18  0.22  0.25  0.29  0.32  0.35  0.39  0.37 ... 
             0.58  0.73  0.96  1.34  2.10  4.39]';
     u = [1:15]';   v = 16 - u;   w = min([u'; v']).';
     xm = [8.24106e-02;  1.13304e+00;  2.34370e+00];
     par = struct('p',8, 'pt',0, 'uvwy',[u v w Bard], 'xm',xm);
     x0 = ones(1,3);   tau0 = 1e-8;   delta0 = 1;
   case 9    % Kowalik and Osborne
     xm = [1.92807e-01;  1.91282e-01;  1.23057e-01; 1.36062e-01];
     y = [0.1957  0.1947  0.1735  0.1600  0.0844  0.0627 ...
          0.0456  0.0342  0.0323  0.0235  0.0246]';
     u = [4.0000  2.0000  1.0000  0.5000  0.2500  0.1670 ...
          0.1250  0.1000  0.0833  0.0714  0.0625]';
     par = struct('p',9, 'pt',0, 'yu',[y u], 'xm',xm);
     x0 = [.25  .39  .415  .39];  tau0 = 1;   delta0 = .1;
   case 10    % Meyer
     xm = [5.6099528e-03;  6.18129930e+03;  3.45222051e+02];
     par = struct('p',10, 'pt',0, 'ty',[45+5*(1:16)' Meyer], 'xm',xm);
     x0 = [.02  4e3  250];   tau0 = 1e-6;   delta0 = 100;
   case 11    % Watson
     t = [1:29]'/29;   B = ones(29,n);
     for  j = 2 : n,  B(:,j) = t.*B(:,j-1); end
     A = B(:,1:n-1) * diag(1:n-1);   xm = [];
     switch  n
       case  6,  xm = [-1.57251e-02;  1.01243e+00; -2.32992e-01
                        1.26043e+00; -1.51373e+00;  9.92996e-01];
       case  9,  xm = [-1.53070e-05;  9.99790e-01;  1.47640e-02
                        1.46342e-01;  1.00082e+00; -2.61773e+00
                        4.10440e+00; -3.14361e+00;  1.05263e+00];
       case 12,  xm = [-6.65800e-09;  1.00000e+00; -5.62881e-04
                        3.47798e-01; -1.56525e-01;  1.05178e+00
                       -3.24414e+00;  7.28248e+00; -1.02647e+01
                        9.06875e+00; -4.53913e+00;  1.01161e+00];
     end
     par = struct('p',11, 'pt',0, 'A',A, 'B',B, 'xm',xm);
     x0 = zeros(1,n);   tau0 = 1e-8;   delta0 = 1;
   case 12    % Box 3-d
     par = struct('p',12, 'pt',0, 't', .1*[1:m]', 'xm',[1; 10; 1]);
     x0 = [0  10  20];     tau0 = 1e-8;   delta0 = 1;
   case 13    % Jennrich and Sampson
     if  m == 10, xm = .257825*[1;1]; else, xm = []; end
     par = struct('p',13, 'pt',0, 't', [1:m]', 'xm',xm);
     x0 = [.3  .4];  tau0 = 1;   delta0 = .05;
   case 14    % Brown and Dennis
     if  m == 20
       xm = [-1.159444e+01;  1.320363e+01; -4.0343945e-01; 2.367787e-01]; 
     else, xm = []; end
     par = struct('p',14, 'pt',0, 't', [1:m]'/5, 'xm',xm);
     x0 = [25  5  -5  -1];   tau0 = 1e-3;   delta0 = .5;
   case 15    % Chebyquad
     y = zeros(m,1);  xm = [];
     if      (n == 8) & (m == 8)
       xm = [4.31528e-02;  1.93091e-01;  2.66329e-01
             5.00000e-01;  5.00000e-01;  7.33671e-01
             8.06909e-01;  9.56847e-01];
     elseif  (n == 8) & (m == 16)
       xm = [7.97431e-02;  1.78355e-01;  3.10120e-01
             4.41205e-01;  5.58795e-01;  6.89880e-01
             8.21645e-01;  9.20257e-01];
     elseif  (n == 9) & (m == 9)
       xm = [4.42053e-02;  1.99491e-01;  2.35619e-01
             4.16047e-01;  5.00000e-01;  5.83953e-01
             7.64381e-01;  8.00509e-01;  9.55795e-01];
     elseif  (n == 9) & (m == 18)
       xm = [1.01891e-01;  1.89483e-01;  2.94482e-01
             3.97099e-01;  5.00000e-01;  6.02901e-01
             7.05518e-01;  8.10517e-01;  8.98109e-01];
     end
     for  i = 2 : 2 : m,  y(i) = -1/(i^2 - 1); end
     par = struct('p',15, 'pt',0, 'y', y, 'xm',xm); 
     x0 = [1:n]/(n+1);  tau0 = 1;   delta0 = 1/(n+1);
   case 16    % Brown almost linear
     A = eye(n-1,n) + ones(n-1,n);
     par = struct('p',16, 'pt',0, 'A', A, 'xm',ones(n,1));
     x0 = .5*ones(1,n);  tau0 = 1;   delta0 = 1;
   case 17    % Osborne 1
     xm = [3.75410e-01;  1.93585e+00; -1.46469e+00
           1.28675e-02;  2.21227e-02];
     par = struct('p',17, 'pt',0, 'ty',Osborne1, 'xm',xm);
     x0 = [.5  1.5  -1  .01  .02];     tau0 = 1e-8;   delta0 = .1;
   case 18    % Exponential fit.  n = 4
     xm = [-4.00003; -4.99996;  4.00025; -4.00025];
     par = struct('p',18, 'pt',0, 'ty',hbnefit, 'xm',xm);
     x0 = [-1 -2 1 -1];     tau0 = 1e-3;   delta0 = 1;
   case 19    % Exponential fit.  n = 2
     xm = [-4.00003; -4.99996];
     par = struct('p',19, 'pt',0, 'ty',hbnefit, 'xm',xm);
     x0 = [-1 -2];     tau0 = 1e-3;   delta0 = 1;
   case 20    % Scaled Meyer
     xm = [2.48178;  6.18135;  3.45224];
     par = struct('p',20, 'pt',0, 'ty',[.45+.05*(1:16)' 1e-3*Meyer], 'xm',xm);
     x0 = [8.85 4 2.5];   tau0 = 1;   delta0 = 1;
   case 21    % Separated Meyer 
     xm = [6.18135e+03;  3.45224e+02];
     par = struct('p',21, 'pt',0, 'ty',[45+5*(1:16)' Meyer], 'xm',xm);
     x0 = [4e3  250];   tau0 = 1;   delta0 = 100;
   case 22    % Exp and squares 
     nn = [1:n].^2;   A = sum(1./nn);
     y = .75;   more = 1;
     while  more
       yo = y;   e = A*exp(-y);   y = y - (y-e)/(1+e);
       more = abs(y-yo) > 1e-12;
     end
     xm = exp(-y)./nn;
     par = struct('p',22, 'pt',1, 'xm',xm);
     x0 = zeros(1,n);   tau0 = 1e-3;   delta0 = 1;
     
   end
   
function  y = Meyer
  y = [34780  28610  23650  19630  16370  13720  11540 ...
        9744   8261   7030   6005   5147   4427   3820 ...
        3307   2872]';
       
function  ty = Osborne1
  y = [0.844  0.908  0.932  0.936  0.925  0.908  0.881 ...
       0.850  0.818  0.784  0.751  0.718  0.685  0.658 ...
       0.628  0.603  0.580  0.558  0.538  0.522  0.506 ...
       0.490  0.478  0.467  0.457  0.448  0.438  0.431 ...
       0.424  0.420  0.414  0.411  0.406]';
  t = 10 * (0 : 32);  ty = [t(:) y];
       
function  ty = hbnefit
  y =[0.090542  0.124569  0.179367  0.195654  0.269707 ...
      0.286027  0.289892  0.317475  0.308191  0.336995 ...
      0.348371  0.321337  0.299423  0.338972  0.304763 ...
      0.288903  0.300820  0.303974  0.283987  0.262078 ...
      0.281593  0.267531  0.218926  0.225572  0.200594 ...
      0.197375  0.182440  0.183892  0.152285  0.174028 ...
      0.150874  0.126220  0.126266  0.106384  0.118923 ...
      0.091868  0.128926  0.119273  0.115997  0.105831 ...
      0.075261  0.068387  0.090823  0.085205  0.067203]';
  t = .02 * [1:45];   ty = [t(:) y];
