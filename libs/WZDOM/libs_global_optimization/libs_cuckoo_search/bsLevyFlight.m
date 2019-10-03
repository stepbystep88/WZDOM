function [nests, fitness] = bsLevyFlight(objFunc, beta, nests, Lb, Ub, a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A subroutine of the Modified Cuckoo Search (MCS) algorithm
% I re-organize this code as a seperate function is because it would be better
% to develop more programs using some common subroutines.
%
% Organized by Bin She (bin.stepbystep@gmail.com)
% Organizing dates: June 2019
% ----------------------------------------------------------------%

    % Levy flights
    [nDim, nNest] = size(nests);
    fitness = zeros(nNest, 1);
    
    for i = 1 : nNest
        nest = nests(:, i);
        step = bsLevy(nDim, 1, beta);
        
%         direction = sign(rand(nDim, 1) - 0.5);
        nest = nest + a .* step;
        nests(:, i) = bsSimpleBounds(nest, Lb, Ub);
        fitness(i) = objFunc(nests(:, i));
    end
end