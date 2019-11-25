function [x, output] = bsPreBuildModelParam(welllog, mode, extraInfo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is to build a model parameter to be estimated.
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
% Input
%
% vp            welllog(:, 2)
% vs            welllog(:, 3)
% rho           welllog(:, 4)
% mode          the flag of how to build the model parameter x
% extraInfo     extra information needed in the building process
% -------------------------------------------------------------------------
% Output
% x             model parameter
% output        some other information
% -------------------------------------------------------------------------
    
    vp = welllog(:, 2);
    vs = welllog(:, 3);
    rho = welllog(:, 4);
    
    Lp = log(vp .* rho);
    Ls = log(vs .* rho);
    Ld = log(rho);
    sampNum = length(Lp);
    output = [];
    
    switch mode
        case 'Lpsd_fit'
            if isfield(extraInfo, 'lsCoef') && ~isempty(extraInfo.lsCoef)
                output = extraInfo;
            else
                output.lsCoef = polyfit(Lp, Ls, 1);
                output.ldCoef = polyfit(Lp, Ld, 1);
            end
            
            output.deltaLs = Ls - ( lsCoef(1)*Lp + lsCoef(2) ); 
            output.deltaLd = Ld - ( ldCoef(1)*Lp + ldCoef(2) );
            x = [Lp; deltaLs; deltaLd];
            
        case 'Lpsd'
            x = [Lp; Ls; Ld];
            
        case 'reflectivity'
            D = stpGen1DDiff(sampNum, 1, 1);
            x = [D*Lp; D*Ls; D*Ld];
            
        otherwise
            validatestring(mode, {'Lpsd_fit', 'Lpsd', 'reflectivity'});
    end
    
end