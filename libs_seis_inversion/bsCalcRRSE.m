function RRSEs = bsCalcRRSE(truemodel, initial, output)
    nTrace = size(truemodel, 2);
    RRSEs = zeros(1, nTrace);
    
    for i = 1 : nTrace
%         residual = ;
%         RMSE = sqrt(mse(output(:, i) - truemodel(:, i)));
        RMSE = norm(output(:, i) - truemodel(:, i), 1);
%         RRMSE = RMSE / (sqrt(mse(initial(:, i)-truemodel(:, i))));
%         RRMSE = RMSE / norm(initial(:, i) - truemodel(:, i), 1);
        RRMSE = RMSE / norm(truemodel(:, i));
        
%         tmp = corrcoef(output(:, i), truemodel(:, i));
%         RRMSE = tmp(1, 2);
        RRSEs(i) = RRMSE;
    end
    
end