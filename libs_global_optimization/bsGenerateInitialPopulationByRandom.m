function population = bsGenerateInitialPopulationByRandom(Lb, Ub, N)
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
    
    population = zeros(nDim, N);

    nDim = length(Lb);
    
    % Random initial solutions
    for i = 1 : N
        population(:, i) = Lb + (Ub - Lb) .* rand(nDim, 1);
    end
    
end