function population = bsSimpleBoundsByChaos(population, Lb, Ub)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code simply projects the population into the range of lower boundary
% Lb and upper boundary Ub
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: May 2019
% -------------------------------------------------------------------------
% Input
%
% population: a set of varibles with size of (the number of dimensions, the
% number of varibles)
%
% Lb: lower boundary
%
% Ub: upper boundary
%
% -------------------------------------------------------------------------
% Output
%
% population: generated population with size (length(Lb), N)
% -------------------------------------------------------------------------

    [nDim, nVar] = size(population);
    
    if ~isempty(Lb)
        % by using repmat, we don't need to write a for loop
        repLb = repmat(Lb, 1, nVar);
        I = population < repLb;
%         population(I) = repLb(I);
    end

    % Apply the upper bounds 
    if ~isempty(Ub)
        repUb = repmat(Ub, 1, nVar);
        I = population > repUb;
%      	population(I) = repUb(I);
    end
    
    seqs = zeros(nDim, nVar);
    for i = 1 : nVar
        seq = bsGenChaos(3, rand(), nDim);
        seqs(:, i) = seq';
    end
    
    population(I) = repLb(I) + (repUb(I) - repLb(I)) .* seqs(I);
end