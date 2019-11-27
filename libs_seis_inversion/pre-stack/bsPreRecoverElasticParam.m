function [vp, vs, rho] = bsPreRecoverElasticParam(x, mode, extraInfo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is to recover elastic parameter from a given model parameter
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------

    sampNum = length(x) / 3;
    nLp = x(1 : sampNum);
    nLs = x( sampNum+1 : sampNum*2);
    nLd = x( sampNum*2+1 : sampNum*3);
    
    estLp = nLp;
    
    switch lower(mode)
        case 'lpsd_fit'
            estLs = extraInfo.lsCoef(1) * nLp + extraInfo.lsCoef(2) + nLs;
            estLd = extraInfo.ldCoef(1) * nLp + extraInfo.ldCoef(2) + nLd;
            
            
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