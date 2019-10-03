function [xOut, fval, exitFlag, output] = bsSeisInv1DByPicewiseSolver(d, G, xInit, Lb, Ub, regParam, parampkgs, options, regFunc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is designed for 1D seismic inversion using picewise
% smooth/linear regularization method
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: May 2019
% -------------------------------------------------------------------------
% Input
% 
% d: observed seismic data, a vector.
% 
% G: forward operator made up by the wavelet information
% 
% xInit: initial guess of model parameters
% 
% Lb: lower boundary of x
% 
% Ub: upper boundary of x
% 
% regParam: regularization parameter. If it is empty, I will start a search
% process to find the optimal regParam.
%
% parampkgs: specail parameter for diffrent methods.
% 
% options: options parameters for 1D seismic inversion. See function
% bsCreateSeisInv1DOptions
% 
% regFunc: regularization function handle. It could be @bsPicewiseSmoothReg
% or @bsPicewiseLinearReg
%
% -------------------------------------------------------------------------
% Output
%
% x          is a column vector; refers to the estimated result.
%
% fval          the objective function value at the last iteration
%
% exitFlag      corresponds to the stopCriteria
% see function bsCheckStopCriteria
%
% output    a struct, in general, it has
% output.iterations: the number of iterations
% output.nfev: the number of function evaluations
% ouput.midResults: the middle results during the iteration process
% ouput.gradient: the gradient at the final estimated parameter
% ouput.options: the options for this inversion process
% output.regParam: the regularization parameters used
% -------------------------------------------------------------------------


    % the struct for calculating residual term 
    mainData.A = G;
    mainData.B = d;
    
    % whether add the low frequency component constraint |x-x0|_2^2
    if options.addLowFreqConstraint
        initObjFunc = {@bsReg1DTKInitModel, [], options.initRegParam};
    else
        initObjFunc = [];
    end
    
    % set default 
    [~, parampkgs] = bsGetFieldsWithDefaults(parampkgs, {
        'diffOrder1', 1; 
        'diffOrder2', 2; 
        'rho', options.initRegParam;
        'nSegments', options.nSegments;
        'D1', [];
        'D2', [];
        });
    
    if isempty(parampkgs.D1) || isempty(parampkgs.D2)
        parampkgs.D1 = bsGen1DDiffOperator(length(xInit), parampkgs.nSegments, parampkgs.diffOrder1);
        parampkgs.D2 = bsGen1DDiffOperator(length(xInit), parampkgs.nSegments, parampkgs.diffOrder2);
        
        % for the first term
        regpkgs.nSegments = parampkgs.nSegments;
        regpkgs.diffOrder = parampkgs.diffOrder1;  
        inputObjFcnPkgs1 = [{options.mainFunc, mainData, 1; @bsReg1DTV, regpkgs, 0};  initObjFunc ];
        
        [~, parampkgs.inputObjFcnPkgs1] = bsCreateFuncPkgs(xInit, inputObjFcnPkgs1);
        
        % for the second term
        regpkgs.diffOrder = parampkgs.diffOrder2;
        switch func2str(regFunc)
            % re-organize the input objective function pakages
            case {'bsPicewiseLinearReg', 'bsPicewiseLinearReg2'}
                inputObjFcnPkgs2 = [{options.mainFunc, mainData, 1; @bsReg1DTV, regpkgs, 0};  initObjFunc ];
            case {'bsPicewiseSmoothReg', 'bsPicewiseSmoothReg2' }
                inputObjFcnPkgs2 = [{options.mainFunc, mainData, 1; @bsReg1DTK, regpkgs, 0};  initObjFunc ];
        end
        [~, parampkgs.inputObjFcnPkgs2] = bsCreateFuncPkgs(xInit, inputObjFcnPkgs2);
    end
    
    if isempty(regParam)
        lambda1 = bsFindBestRegParameter(options, parampkgs.inputObjFcnPkgs1, xInit, Lb, Ub);
        lambda2 = bsFindBestRegParameter(options, parampkgs.inputObjFcnPkgs2, xInit, Lb, Ub);
        
    else
        lambda1 = regParam.lambda1;
        lambda2 = regParam.lambda2;
    end

    % set the regularization parameters
    parampkgs.lambda1 = lambda1;
    parampkgs.lambda2 = lambda2;
    parampkgs.inputObjFcnPkgs1{2, 3} = lambda1;
    parampkgs.inputObjFcnPkgs2{2, 3} = lambda2;
    parampkgs.rho = options.initRegParam;
   
    % set the main function made up by the residual term and the low 
    % frequency componet constraint (if it requires)
    inputMainFcnpkgs = bsCreateFuncPkgs(xInit, [{options.mainFunc, mainData, 1}; initObjFunc]);
    mainFcnpkgs = @(x)(bsSumFuncsWithBound(x, options.GBOptions.bounds, inputMainFcnpkgs));
            
    % prepare initial model
    GBOptions = options.GBOptions;
    GBOptions.maxIter = 100;
    GBOptions.display = 'off';
    GBOptions.plotFcn = [];
    GBOptions.isSaveMiddleRes = false;
    parampkgs.x1 = bsGBSolveByOptions(parampkgs.inputObjFcnPkgs1, xInit, [], [], GBOptions); 
    parampkgs.x2 = bsGBSolveByOptions(parampkgs.inputObjFcnPkgs2, xInit, [], [], GBOptions); 
    
    % call regFunc to create packages and function handles to excute in
    % bsMultiParamSolver
    [pkgs, fcnUpdateParams, fcnGetInitialValue] = regFunc(xInit, Lb, Ub, mainFcnpkgs, parampkgs);
    
    % call MultiParamSolver to solve the problem
    [xOut, fval, exitFlag, output] ...
        = bsMultiParamSolver(fcnUpdateParams, fcnGetInitialValue, pkgs, options.GBOptions, options.maxIter, options.innerIter);
    
    % save the regularization parameter
    regParam.lambda1 = lambda1;
    regParam.lambda2 = lambda2;
    
    output.regParam = regParam;
    output.parampkgs = parampkgs;
end
