function [welllog] = bsFiltWelllog(welllog, filtCoef)

    for i = 2 : 4
        welllog(:, i) = bsButtLowPassFilter(welllog(:, i), filtCoef);
    end
%     depth = welllog(:, indexDepth);
%     vp = bsButtLowPassFilter(welllog(:, indexVp), filtCoef);
%     vs = bsButtLowPassFilter(welllog(:, indexVs), filtCoef);
%     rho = bsButtLowPassFilter(welllog(:, indexRho), filtCoef);
end