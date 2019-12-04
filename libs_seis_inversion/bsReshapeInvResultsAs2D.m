function invResults = bsReshapeInvResultsAs2D(invResults)
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
    
    for i = 1 : length(invResults)
        data = invResults{i}.data;
        
        
        if ~iscell(data)
            fprintf('Reshaping %s data of %s...\n', invResults{i}.type, invResults{i}.name);
            invResults{i}.data = bsReshapeData(invResults{i}.data);
        else
            for j = 1 : length(data)
                fprintf('Reshaping %s data of %s...\n', invResults{i}.type{j}, invResults{i}.name);
                invResults{i}.data{j} = bsReshapeData(invResults{i}.data{j});
            end
        end
        
    end
end

function volume = bsReshapeData(data)
    [sampNum, nCrossline, nInline] = size(data);
    volume = reshape(data, sampNum, nCrossline*nInline);
end