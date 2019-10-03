function newNests = bsEmptyNestsWithAlpha(nests, Lb, Ub, pa, alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A subroutine of the Cuckoo Search (CS) algorithm
% To Replace some nests by constructing new solutions/nests
%
% I re-organize this code as a seperate function is because it would be better
% to develop more programs using some common subroutines.
%
% Organized by Bin She (bin.stepbystep@gmail.com)
% Organizing dates: May 2019
% 
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

    [nDim, nNest] = size(nests);
    
    alphas = normrnd(alpha.value, 0.1, 1, nNest);
%     alphas = rand(nNest, 1);
    alphas(alphas > alpha.max) = alpha.max;
    alphas(alphas < alpha.min) = alpha.min;
    alphas = repmat(alphas, nDim, 1);
    
    % Discovered or not -- a status vector
    K = rand(size(nests)) <= pa;
    % In the real world, if a cuckoo'nest egg is very similar to a host'nest eggs, then 
    % this cuckoo'nest egg is less likely to be discovered, thus the fitness should 
    % be related to the difference in solutions.  Therefore, it is a good idea 
    % to do a random walk in a biased way with some random step sizes.  
    
    %% New solution by biased/selective random walks
    stepsize = alphas .* (nests(:, randperm(nNest)) - nests(:, randperm(nNest)));
    newNests = nests + stepsize .* K;
    
    % call bsSimpleBounds to project the new nests into the reasonable range
    newNests = bsSimpleBounds(newNests, Lb, Ub);
    
end