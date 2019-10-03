function [population] = bsGenerateInitialPopulationByLHS(Lb, Ub, N)
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
%     nDim = length(Lb);
    
%     population = zeros(nDim, N);
%     diff = Ub - Lb;

%     lhsSamples = lhsdesign(nDim, N);
% 
%     repLb = repmat(Lb, 1, N);
%     repDiff = repmat(diff, 1, N);
    
%     population = repLb + repDiff .* lhsSamples;
    population = bsLHS(Lb, Ub, N);
    
end

