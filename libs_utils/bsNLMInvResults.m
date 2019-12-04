function invResults = bsNLMInvResults(invResults, refData, varargin)
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
    
    weightInfo = [];
    
    for i = 1 : length(invResults)
        data = invResults{i}.data;
        
        if ~iscell(data)
            fprintf('Smoothing %s data of %s by using NLM filter.\n', invResults{i}.type, invResults{i}.name);
            [invResults{i}.data, weightInfo] = bsNLMByRef(data, refData, ...
                'weightInfo', weightInfo, varargin{:});
        else
            for j = 1 : length(data)
                fprintf('Smoothing %s data of %s by using NLM filter.\n', invResults{i}.type{j}, invResults{i}.name);
                [invResults{i}.data{j}, weightInfo] = bsNLMByRef(data{j}, refData, ...
                    'weightInfo', weightInfo, varargin{:});
            end
        end
        
    end
end