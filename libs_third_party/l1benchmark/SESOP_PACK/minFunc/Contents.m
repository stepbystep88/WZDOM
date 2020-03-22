% MINFUNC: L-BFGS and many other methods, Mark Schmidt (2006)
%
% Files
%   ArmijoBacktrack             - Backtracking linesearch to satisfy Armijo condition
%   autoGrad                    - [f,g] = autoGrad(x,useComplex,funObj,varargin)
%   autoHess                    - Numerically compute Hessian of objective function from gradient values
%   autoHv                      - Numerically compute Hessian-vector product H*v of funObj(x,varargin{:})
%   autoTensor                  - Numerically compute Tensor of 3rd-derivatives of objective function from Hessian values
%   conjGrad                    - Linear Conjugate Gradient, where optionally we use
%   example_minFunc             - Runs various limited-memory solvers on 2D rosenbrock function for 25
%   lbfgs                       - BFGS Search Direction
%   lbfgsUpdate                 - 
%   mchol                       - Compute a modified LDL factorization of A
%   mcholinc                    - Computes Cholesky of H+tau*I, for suitably large tau that matrix is pd
%   minFunc                     - minFunc(funObj,x0,options,varargin)
%   minFunc_processInputOptions - 
%   minFuncBC                   - minFunc(funObj,x0,LB,UB,options,varargin)
%   minFuncEC                   - minFunc(funObj,x0,Aeq,Beq,options,varargin)
%   minFuncIC                   - minFunc(funObj,x0,A,B,options,varargin)
%   polyinterp                  - Minimum of interpolating polynomial based on function and derivative
%   precondDiag                 - 
%   precondTriu                 - 
%   precondTriuDiag             - 
%   rosenbrock                  - Calculate objective f
%   taylorModel                 - 
%   WolfeLineSearch             - Bracketing Line Search to Satisfy Wolfe Conditions
