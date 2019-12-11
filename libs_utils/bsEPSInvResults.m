function invResults = bsEPSInvResults(invResults, varargin)
%% smooth data by using the EPS algorithm
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Dec 2019
% 
% Input
% options.windowSize        the size of the sliding 2D window

% see bsEPSSmooth.m function
% -------------------------------------------------------------------------

    if isempty(invResults)
        return;
    end
    
    for i = 1 : length(invResults)
        data = invResults{i}.data;
        
        if ~iscell(data)
            fprintf('Smoothing %s data of %s by using EPS filter.\n', invResults{i}.type, invResults{i}.name);
            invResults{i}.data = bsEPSmooth(data, varargin{:});
        else
            for j = 1 : length(data)
                fprintf('Smoothing %s data of %s by using EPS filter.\n', invResults{i}.type{j}, invResults{i}.name);
                invResults{i}.data{j} = bsEPSmooth(data{j}, varargin{:});
            end
        end
        
    end
end