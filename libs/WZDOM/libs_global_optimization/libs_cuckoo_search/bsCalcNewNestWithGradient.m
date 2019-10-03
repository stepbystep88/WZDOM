function [bestFVal, bestNest, newNest, fitness, gradient, K] = bsCalcNewNestWithGradient(fobj, newNest)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A subroutine of the Cuckoo Search (CS) algorithm
% To find the current bestNest nests with gradient information 
%
% I re-organize this code as a seperate function is because it would be better
% to develop more programs using some common subroutines.
%
% Organized by Bin She (bin.stepbystep@gmail.com)
% Organizing dates: May 2019
% -----------------------------------------------------------------
% Cuckoo Search (CS) algorithm by Xin-She Yang and Suash Deb      %
% Programmed by Xin-She Yang at Cambridge University              %
% Programming dates: Nov 2008 to June 2009                        %
% Last revised: Dec  2009   (simplified version for demo only)    %
% -----------------------------------------------------------------
% Papers -- Citation Details:
% 1) X.-S. Yang, S. Deb, Cuckoo search via Levy flights,
% in: Proc. of World Congress on Nature & Biologically Inspired
% Computing (NaBIC 2009), December 2009, India,
% IEEE Publications, USA,  pp. 210-214 (2009).
% http://arxiv.org/PS_cache/arxiv/pdf/1003/1003.1594v1.pdf 
% 2) X.-S. Yang, S. Deb, Engineering optimization by cuckoo search,
% Int. J. Mathematical Modelling and Numerical Optimisation, 
% Vol. 1, No. 4, 330-343 (2010). 
% http://arxiv.org/PS_cache/arxiv/pdf/1005/1005.2908v2.pdf
% ----------------------------------------------------------------%

    gradient = zeros(size(newNest));
    fitness = zeros(size(newNest, 2), 1);
    
    % Evaluating all new solutions
    for j = 1 : size(newNest, 2)
        [fNew, gNew] = fobj(newNest(:, j));
        fitness(j) = fNew;
        gradient(:, j) = gNew;
    end
    
    % Find the current bestNest
    [bestFVal, K] = min(fitness) ;
    bestNest = newNest(:, K);
end
