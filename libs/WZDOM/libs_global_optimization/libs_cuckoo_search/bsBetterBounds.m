function new_population = bsBetterBounds(trial_population, Lb, Ub, old_population)
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

    [nDim, nVar] = size(trial_population);
    new_population = trial_population;
    
    % Apply the lower bounds 
    if ~isempty(Lb)
        % by using repmat, we don't need to write a for loop
        repLb = repmat(Lb, 1, nVar);
        
        I = trial_population < repLb;
        % it means new_population = 2*Lb - trial_population for the elements of old_population
        % that are smaller than lower boundary
        new_population(I) = 2 * repLb(I) - trial_population(I);   
    end

    % Apply the upper bounds 
    if ~isempty(Ub)
        repUb = repmat(Ub, 1, nVar);
        
        I = trial_population > repUb;
        % it means new_population = 2*Ub - new_population for the elements of trial_population
        % that are greater than upper boundary
        new_population(I) = 2 * repUb(I) - trial_population(I);
    end

    if ~isempty(old_population)
        I = new_population < repLb;
        % we then set the new_population = 0.5 * (lower + old_population),
        % note that the old_population satisfies the boudary constraints
        new_population(I) = (repLb(I) + old_population(I)) * 0.5;    
    end

    if ~isempty(old_population)
        I = new_population > repUb;
        % we then set the new_population = 0.5 * (lower + old_population),
        % note that the old_population satisfies the boudary constraints
        new_population(I) = (repUb(I) + old_population(I)) * 0.5;
    end
end