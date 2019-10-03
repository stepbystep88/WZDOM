function [population] = bsGenerateInitialPopulationByChaos(Lb, Ub, N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is designed for genearating latin hypercube sampling 
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: May 2019
% -------------------------------------------------------------------------
% Input
%
% Lb: lower boundary
%
% Ub: upper boundary
% 
% N: the number of population
%
% -------------------------------------------------------------------------
% Output
%
% population: generated population with size (length(Lb), N)
% -------------------------------------------------------------------------


    % the number of dimensions
    nDim = length(Lb);
    
%     population = zeros(nDim, N);
    diff = Ub - Lb;

    
%     lhsSamples = lhsdesign(nDim, N, 'smooth', 'off', 'criterion', 'maximin', 'iterations', 100);

    population = zeros(nDim, N);
    
    for i = 1 : N
        seq = bsGenChaos(3, rand(), nDim);
        population(:, i) = Lb + diff .* seq';
    end
    
%     population = repLb + repDiff .* lhsSamples;
    
end

