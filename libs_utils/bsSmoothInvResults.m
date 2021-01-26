function invResults = bsSmoothInvResults(invResults, refData, fcn, varargin)
%% smooth data by using the NLM algorithm, the similarites are referenced from refData
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Dec 2019
% 
% Input
% options.p                 parameter for inverse distance weight function,
%                           its reasonable range is [0.5 3]
% options.nPointsUsed       the most number of points used to calculate the
%                           weight information
% options.stride            step size of sliding window
% options.searchOffset      indicating the search range
% options.windowSize        the size of the sliding 2D window

% see bsNLMByRef.m function
% -------------------------------------------------------------------------

    if isempty(invResults)
        return;
    end
    
    rangeInline = [min(invResults{1}.inIds), max(invResults{1}.inIds)];
    rangeCrossline = [min(invResults{1}.crossIds), max(invResults{1}.crossIds)];
    
    nCrossline = rangeCrossline(2) - rangeCrossline(1) + 1;
    nInline = rangeInline(2) - rangeInline(1) + 1;
    
    if isempty(fcn)
        fcn = @bsSmoothByGST2D;
    end
    
    weightInfo = [];
    
    for i = 1 : length(invResults)
        data = invResults{i}.data;
        
        if ~iscell(data)
            fprintf('Smoothing %s data of %s by using %s filter.\n', invResults{i}.type, invResults{i}.name, char(fcn));
            
            if strcmpi(char(fcn), 'bsNLMByRef')
                [invResults{i}.data, weightInfo] = bsNLMByRef(data, refData, ...
                    'weightInfo', weightInfo, varargin{:});
            else
                if length(invResults{i}.inIds) > 10000
                    % 表明这是一个3D数据
                    data = bsReshapeDataAs3D(data, nInline, nCrossline);
                    tmp = fcn(data, refData, varargin{:});
                    invResults{i}.data = bsReshapeDataAs2D(tmp);
                else
                    [invResults{i}.data] = fcn(data, refData, varargin{:});
                end
            end
            
        else
            for j = 1 : length(data)
                fprintf('Smoothing %s data of %s by using %s filter.\n', invResults{i}.type{j}, invResults{i}.name, char(fcn));
                if strcmpi(char(fcn), 'bsNLMByRef')
                    [invResults{i}.data{j}, weightInfo] = bsNLMByRef(data{j}, refData, ...
                        'weightInfo', weightInfo, varargin{:});
                else
                    if length(invResults{i}.inIds) > 20000
                        % 表明这是一个3D数据
                        tmp1 = bsReshapeDataAs3D(data{j}, nInline, nCrossline);
                        if ~isempty(refData) && length(size(refData)) < 3
                            refData = bsReshapeDataAs3D(refData, nInline, nCrossline);
                        end
                        
                        tmp = fcn(tmp1, refData, varargin{:});
                        invResults{i}.data{j} = bsReshapeDataAs2D(tmp);
                    else
                        [invResults{i}.data{j}] = fcn(data{j}, refData, varargin{:});
                    end
                
                    
                end
            end
        end
        
    end
end