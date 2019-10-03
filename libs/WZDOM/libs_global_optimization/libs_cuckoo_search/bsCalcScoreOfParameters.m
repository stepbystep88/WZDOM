function algParams = bsCalcScoreOfParameters(algParams, fitness, newFitNess, flag, nHistory)

    delatFit = fitness - newFitNess;
    increasedFit = sum(delatFit(delatFit>0)) / length(fitness);
%     increasedFit = sum(delatFit>0);

%     if increasedFit > 0
        switch flag
            case 1
                algParams.alpha1 = bsUpdateParameter(algParams.alpha1, increasedFit, nHistory);
                algParams.alpha2 = bsUpdateParameter(algParams.alpha2, increasedFit, nHistory);
                algParams.Fi = bsUpdateParameter(algParams.Fi, increasedFit, nHistory);
                algParams.pe = bsUpdateParameter(algParams.pe, increasedFit, nHistory);
                algParams.pb = bsUpdateParameter(algParams.pb, increasedFit, nHistory);
                algParams.cr = bsUpdateParameter(algParams.cr, increasedFit, nHistory);
                algParams.lambda = bsUpdateParameter(algParams.lambda, increasedFit, nHistory);
                
%                 algParams.alpha.history(iHistory, 2) = increasedFit;
%                 algParams.Fi.history(iHistory, 2) = increasedFit;
%                 algParams.pe.history(iHistory, 2) = increasedFit;
%                 algParams.pb.history(iHistory, 2) = increasedFit;
%                 algParams.cr.history(iHistory, 2) = increasedFit;
%                 algParams.lambda.history(iHistory, 2) = increasedFit;
            case 2
%                 algParams.alpha.pa(iHistory, 2) = increasedFit;
                algParams.pa = bsUpdateParameter(algParams.pa, increasedFit, nHistory);
                algParams.pd = bsUpdateParameter(algParams.pd, increasedFit, nHistory);
        end
%     end
    
end

function t = bsUpdateParameter(t, increasedFit, nHistory)
    t.index = t.index + 1;
    t.nTrain = t.nTrain + 1;
    
    if (t.index > nHistory)
        t.index = 1;
    end
    
    if (t.nTrain > nHistory)
        t.nTrain = nHistory;
    end
    
    t.store(t.index, 2) = increasedFit;
    t.store(t.index, 1) = t.value;
end