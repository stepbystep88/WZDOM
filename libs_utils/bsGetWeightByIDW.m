function [weights, indexies] = bsGetWeightByIDW(xs, ys, xs_know, ys_know, options)
%% calculates the weight information by using inverse distance weight (IDW) algorithm 
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Dec 2019
% 
% Input
% options.p                 its reasonable range is [0.5 3]
% options.nPointsUsed       the most number of points used to calculate the
% weight information
%
% Output
% weights                   the calculated weights of each point [xs(i), ys(i)]
% indexies                  the index of known points used for interpolation
% for i-th data [xs(i) ys(i)], its value will be sum(know_data(:, index) .* weights(:, i))
% -------------------------------------------------------------------------

    if nargin < 5
        options = bsSetFields([], {'nPointsUsed', 4; 'p', 2; });
    end
   
    nTrace = length(xs);
    
    fprintf('Calculating the weight information...\n');
    
     % at most use 4 wells for interpolation
    nPointsUsed = min(options.nPointsUsed, length(xs_know));
    indexies = zeros(nTrace, nPointsUsed);
    weights = zeros(nPointsUsed, nTrace);
    
    for i = 1 : nTrace
        if mod(i, 10000) == 0
            fprintf('Calculating the weight information of trace %d/%d...\n', i, nTrace);
        end
        
        D = sqrt((xs(i)-xs_know).^2 + (ys(i)-ys_know).^2);
        [~, index] = bsMinK(D, nPointsUsed); 
        indexies(i, :) = index';
        
        [minD, minIndex] = min(D);
        if minD == 0
            weights(1, i) = 1;
        else
            weights(:, i) = D(index).^(-options.p);
        end
        
        weights(:, i) = weights(:, i) / sum(weights(:, i));
        
    end
    
end