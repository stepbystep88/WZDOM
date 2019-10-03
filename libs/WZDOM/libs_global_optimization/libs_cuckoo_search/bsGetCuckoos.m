function nests = bsGetCuckoos(beta, nests, bestNest, Lb, Ub, alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A subroutine of the Cuckoo Search (CS) algorithm
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

    % Levy flights
    [nDim, nNest] = size(nests);

    %% Levy flights by Mantegna'nest algorithm
%     steps = bsLevy(size(nests, 1), nNest, beta);
    
    steps = bsLevy(nDim, nNest, beta);
    for j = 1 : nNest
        nest = nests(:, j);
        step = steps(:, j);
        
        % In the next equation, the difference factor (nest-bestNest) means that 
        % when the solution is the bestNest solution, it remains unchanged.     
        stepsize = alpha .* step .* (nest - bestNest);
        % Here the factor 0.01 comes from the fact that L/100 should the typical
        % step size of walks/flights where L is the typical lenghtscale; 
        % otherwise, Levy flights may become too aggresive/efficient, 
        % which makes new solutions (even) jump out side of the design domain 
        % (and thus wasting evaluations).
        % Now the actual random walks or flights
%         nest = nest + stepsize .* randn(size(nest));
        nest = nest + stepsize;
        % Apply simple bounds/limits
        nests(:, j) = bsSimpleBounds(nest, Lb, Ub);
    end
end