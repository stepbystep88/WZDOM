function [vp, vs, rho] = bsPreRecoverElasticParam(x, mode, extraInfo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is to recover elastic parameter from a given model parameter
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------

    sampNum = size(x, 1) / 3;
    nLp = x(1 : sampNum, :);
    nLs = x( sampNum+1 : sampNum*2, :);
    nLd = x( sampNum*2+1 : sampNum*3, :);
    
    if iscell(extraInfo)
        tmp = cell2mat(extraInfo);
        ts = reshape([tmp.lsCoef], 2, []);
        td = reshape([tmp.ldCoef], 2, []);
        ls1 = repmat(ts(1, :), sampNum, 1);
        ls2 = repmat(ts(2, :), sampNum, 1);
        ld1 = repmat(td(1, :), sampNum, 1);
        ld2 = repmat(td(2, :), sampNum, 1);
    else
        ls1 = extraInfo.lsCoef(1);
        ls2 = extraInfo.lsCoef(2);
        ld1 = extraInfo.ldCoef(1);
        ld2 = extraInfo.ldCoef(2);
    end
    
    estLp = nLp;
    
    switch lower(mode)
        case 'lpsd_fit'
            estLs = ls1 .* nLp + ls2 + nLs;
            estLd = ld1 .* nLp + ld2 + nLd;
            
            
        case 'lpsd'
            estLs = nLs;
            estLd = nLd;
            
        otherwise
            validatestring(mode, {'Lpsd_fit', 'Lpsd'});
    end
    
    vp = exp(estLp - estLd);
    vs = exp(estLs - estLd);
    rho = exp(estLd);
    
end