function [options] = bsCreateGInvParam(flag, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a inverse struct guiding the inversion process
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Dec 2019
% -------------------------------------------------------------------------
    p = inputParser;
   
    GSegyInfo = bsCreateGSegyInfo();
    p = bsAddCommonParam(p, GSegyInfo);
    
    switch lower(flag)
        case {'prestack', 'pre-stack'}
            % pre-stack data information, see bsPrePrepareModel.m for
            % details
            addParameter(p, 'preSeisData', struct(...
                'mode', 'offset_one_file', ...
                'separates', [], ...
                'fileName', '', ...
                'GSegyInfo', GSegyInfo ...
            ));
        
            % angle data for angle gather seismic data
            addParameter(p, 'angleData', []);
            
            % the maximum of angles
            addParameter(p, 'maxAngle', 40);
            
            % which framework is used for prestack inversion
            % see bsPreBuildModelParam.m for details
            addParameter(p, 'mode', 'lpsd_fit');
            
            % the fitting coeficient
            % see Simultaneous inversion of pre-stack seismic data, Hampson and Russel (2005) for details
            addParameter(p, 'lsdCoef', []);
            
            % whether to set the delat component (deltaLp, deltaLd) as 0
            % see Simultaneous inversion of pre-stack seismic data, Hampson and Russel (2005) for details
            addParameter(p, 'isInitDeltaZero', 0);
            
            % the number of super traces of each gather, the traces with
            % big offset will be discared (why we use old and new)
            addParameter(p, 'oldSuperTrNum', 15);
            addParameter(p, 'newSuperTrNum', 13);
    
            % the number of angle traces
            addParameter(p, 'angleTrNum', 10);
            
        case {'poststack', 'post-stack'}
            
        otherwise
            validatestring(lower(flag), {'prestack', 'pre-stack', 'poststack', 'post-stack'});
    end
    
    p.parse(varargin{:});  
    options = p.Results;
    options.flag = lower(flag);
end

function [p] = bsAddCommonParam(p, GSegyInfo)
    
    % how many sample points above a horizon is inversed (include the current point)
    addParameter(p, 'upNum', 0);
    
    % how many sample points below a horizon is inversed
    addParameter(p, 'downNum', 0);
    
    % the post-stack data information, this one is not neccesaray for 
    % prestack inversion
    addParameter(p, 'postSeisData', struct(...
        'segyInfo', GSegyInfo, ...
        'fileName', '', ...
        'shiftFileName', '' ...
    ));
    
    % the sample interval of the current work area, in ms
    addParameter(p, 'dt', []);
    
    % whether to read the intermediate models from local 
    addParameter(p, 'isReadMode', 0);
    
    % whether to save the intermediate models
    addParameter(p, 'isSaveMode', 0);
    
    % the path to save the middle or final results
    addParameter(p, 'modelSavePath', './');
    
    % whether to normalize the data d, Gm (deviding by norm(d))
    addParameter(p, 'isNormal', 1);
    
    % which horizon is used 
    addParameter(p, 'usedTimeLineId', 1);
    
    % whether to use parallel function
    addParameter(p, 'isParallel', 1);
    
    % how many number of workers are used to perform a parallel program
    addParameter(p, 'numWorkers', bsGetMaxNumWorkers());
    
    % whether to print the exact progress inforamtion when perform a parallel
    % program. If you want to see the progress information, set it as 1. But it should
    % be noted that this setting will decrease the efficiency a little bit and 
    % call I/O many times.
    addParameter(p, 'isPrintBySavingFile', 0);
    
    % the main frequency of wavelet
    addParameter(p, 'waveletFreq', 45);
    
    % wavelet information
    addParameter(p, 'wavelet', []);
    
    % the coeficient of a low-pass filter used for processing seismic data
    addParameter(p, 'seismicFiltCoef', 1);
    
    % the index of each attribute save in welllog data
    addParameter(p, 'indexInWellData', struct(...
        'time', 0, ...
        'ip', 0, ...
        'depth', 1, ...
        'vp', 2, ...
        'vs', 3, ...
        'rho', 4 ...
    ));

    % match depth information and time information
    addParameter(p, 'depth2time', struct(...
        'showCompareNum', 10, ...
        'isShowCompare', 1, ...
        'searchOffsetNum', 20, ...
        'expandNum', 100, ...
        'saveOffsetNum', 20 ...
    ));

    % the function handle of finding regularization parameter
    addParameter(p, 'searchRegParamFcn', @bsBestParameterByLCurve);
    
    % the information of initial model
    info = struct('segyInfo', GSegyInfo, 'fileName', '');
    addParameter(p, 'initModel', struct(...
        'mode', 'filter_from_true_log', ... % see bsPrePrepareModel for details
        'filtCoef', 0.1, ...
        'vp', info, ...
        'vs', info, ...
        'rho', info, ...
        'ip', info ...
    ));
    
    addParameter(p, 'bound', struct(...
        'mode', 'off', ...
        ... % this is used when mode is 'based_on_init'
        'offset_init', struct('vp', 1000, 'vs', 700, 'rho', 0.5, 'ip', 2000), ...
        ... % the following two entries are used when mode is 'fixed'
        'Lb', [], ...
        'Ub', [] ...
    ));
    
    % add error model
    addParameter(p, 'errorModel', struct(...
        'segyInfo', GSegyInfo, ...
        'fileName', '', ...
        'isUse', 0 ...
    ));

    %% set the options of inversion process
    seisInvOptions = bsCreateSeisInv1DOptions();
    seisInvOptions.GBOptions.optAlgHandle = @bsOptQCG;          % using Quasi-Newton conjugate gradient optimizer
    seisInvOptions.GBOptions.display = 'off';
    seisInvOptions.GBOptions.optAlgParam.updateFlag = 'FR';     % PR formulation, could be also FR, HS, DY
    seisInvOptions.GBOptions.isSaveMiddleRes = false;           % Whether save the middle results
    seisInvOptions.GBOptions.stepTolerance = 1e-8;
    seisInvOptions.GBOptions.optimalGradientTolerance = 1e-10;
    seisInvOptions.GBOptions.functionTolerance = 1e-10;
    seisInvOptions.GBOptions.optimalityTolerance = 1e-10;

    seisInvOptions.addLowFreqConstraint = true;
    seisInvOptions.initRegParam = 0.01;
    seisInvOptions.searchRegParamFcn = [];

    addParameter(p, 'seisInvOptions', seisInvOptions);
    
end