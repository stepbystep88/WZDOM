function newNests = bsGetCuckoosByTwoStrategies(nests, algParams, pbest, Lb, Ub)
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
    
    alphas = normrnd(algParams.alpha1.value, 0.1, nNest, 1);
%     alphas = rand(nNest, 1);
    alphas(alphas > algParams.alpha1.max) = algParams.alpha1.max;
    alphas(alphas < algParams.alpha1.min) = algParams.alpha1.min;

    Fis = normrnd(algParams.Fi.value, 0.1, nNest, 1);
    Fis(Fis > algParams.Fi.max) = algParams.Fi.max;
    Fis(Fis < algParams.Fi.min) = algParams.Fi.min;
    
    newNests = nests;
    
%     pbestNest = nests(:, 1);
    steps = bsLevy(nDim, nNest, algParams.lambda.value);
    
    for i = 1 : nNest
        nest = nests(:, i);
        
        % levy flight strategy
%         step = bsLevy(nDim, 1, algParams.lambda.value);
% %             stepsize = alphas(i) .* step;
% %         pbestNest = nests(:, 1);
%         stepsize = alphas(i) .* step .* (nest - pbestNest);
            
%         choose method
        if rand() < algParams.pb.value
            if rand() < algParams.pe.value
                pbestNest = nests(:, randi(pbest));
            else
                pbestNest = nests(:, randi(nNest));
            end
            
            % levy flight strategy
            
%             stepsize = alphas(i) .* step;
            stepsize = alphas(i) .* steps(:, i) .* (nest - pbestNest);
            
        else
            % randomly choose one pbest from nests (nests has been ordered by fitness)
            if rand() < algParams.pe.value
                pbestNest = nests(:, randi(pbest));
            else
                pbestNest = nests(:, randi(nNest));
            end
            
            r1 = randi(nNest);
            r2 = randi(nNest);
            
            stepsize = Fis(i) * (pbestNest - nest + nests(:, r1) - nests(:, r2));
        end
        
        
        newNests(:, i) = nest + stepsize;
    end
    
    % boundary 
    newNests = bsBetterBounds(newNests, Lb, Ub, nests);
%     newNests = bsSimpleBounds(newNests, Lb, Ub);
    
%     random copy the original nests
    K1 = rand(size(nests)) <= algParams.cr.value;
    K2 = repmat((1:nDim)', 1, nNest) == randi(nDim, nDim, nNest);
    K = K1 | K2;
%     
    newNests = newNests .* K + nests .* (~K);
    
end