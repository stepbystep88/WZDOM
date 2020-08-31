function bsShowIterationProcess(GShowProfileParam, invVals, outputs, wellLog)
    
    nItems = length(invVals);
    model = invVals{1}.model;
    
    
    
    colorTbl = bsGetColormap('separate');
    names = cell(1, nItems);
    
    figure;
    for i = 1 : nItems
        midResults = outputs{i}.midResults;
        nIter = size(midResults.f, 2);
        names{i} = invVals{i}.name;
        
        rrsem = bsCalclRRSEM(midResults.x, model.trueLog);
        
        plot(1:nIter, rrsem, '-', 'color', colorTbl{i}, 'linewidth', 2); hold on;
    end
    
    xlabel('Iteration');
    ylabel('RMSE of model');
    legend(names, 'fontweight', 'bold');
    title(wellLog.name, 'fontweight', 'bold');
    set(gca, 'xlim', [1, nIter]);
    
    bsSetDefaultPlotSet(GShowProfileParam.plotParam);
    
    figure;
    for i = 1 : nItems
        midResults = outputs{i}.midResults;
        nIter = size(midResults.f, 2);
        names{i} = invVals{i}.name;
        plot(1:nIter, midResults.f, '-', 'color', colorTbl{i}, 'linewidth', 2); hold on;
    end
    
    xlabel('Iteration');
    ylabel('Objective function value');
    legend(names, 'fontweight', 'bold');
    title(wellLog.name, 'fontweight', 'bold');
    set(gca, 'xlim', [1, nIter]);
    bsSetDefaultPlotSet(GShowProfileParam.plotParam);
    
end

function rrsem = bsCalclRRSEM(xs, trueLog)
    ips = exp(xs);
    nIter = size(xs, 2);
    rrsem = zeros(1, nIter);
    
    for i = 1 : nIter
        ip = ips(:, i);
        rrsem(i) = sqrt(mse((ip - trueLog)));
    end
    
end