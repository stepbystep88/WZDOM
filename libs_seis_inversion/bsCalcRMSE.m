function RRSEs = bsCalcRMSE(x1, x2)
    nTrace = size(x1, 2);
    RRSEs = zeros(1, nTrace);
    
    parfor i = 1 : nTrace
        RRSEs(i) = sqrt(mse(x1(:, i)-x2(:, i)));
    end
    
end