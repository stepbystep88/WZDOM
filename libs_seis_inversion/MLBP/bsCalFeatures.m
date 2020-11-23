function features = bsCalFeatures(mainFunc, regFunc, xInit, mainData, regData)
    f1 = mainFunc(xInit, mainData);
    f2 = regFunc(xInit, regData, 0);

    pN = mainData.B;
    residual = mainData.A * xInit - mainData.B;

    baseD = [norm(pN), mean(pN), var(pN)];
    kurtosisD = kurtosis(pN);
    skewnessD = skewness(pN);
%     gevFitD = gevfit(pN);
%                                 normFitD = normfit(pN);


    baseE = [norm(residual), mean(residual), var(residual)];
    kurtosisE = skewness(residual);
    skewnessE = kurtosis(residual);
%     gevFitE = gevfit(residual);
%                                 normFitE = normfit(residual);

%     features = [
%         baseD, kurtosisD, skewnessD, gevFitD, ...
%         baseE, kurtosisE, skewnessE, gevFitE, ...
%         f1, f2, norm(xInit)];

    features = [
        baseD, kurtosisD, skewnessD, ...
        baseE, kurtosisE, skewnessE, ...
        f1, f2, norm(xInit)];

end