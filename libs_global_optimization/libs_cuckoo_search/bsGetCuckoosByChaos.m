function nests = bsGetCuckoosByChaos(beta, nests, bestNest, Lb, Ub)
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
    
%     a = zeros(nNest, 1);
    
%     seq = bsGenChaos(3, rand(), nNest);
% 
    
    
%     for i = 1 : nNest
%         if seq(i) > 1
%             a(i) = 1;
%         elseif seq(i) < 0.05
%             a(i) = 0.05;
%         else
%             a(i) = seq(i);
%         end
%     end
%     a = rand(nNest, 1);
    a = rand(nNest, 1);
    a(a < 0.05) = 0.05;
    steps =  bsLevy(nDim, nNest, beta);
    
    for i = 1 : nNest
        nest = nests(:, i);
        step = steps(:, i);
        
        stepsize = a(i) .* step .* (nest - bestNest);
        
        nest = nest + stepsize;
        nests(:, i) = bsBetterBounds(nest, Lb, Ub, nests);
    end
end