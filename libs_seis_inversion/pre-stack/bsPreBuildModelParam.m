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
    
    switch lower(mode)
        case 'lpsd_fit'
            if isfield(extraInfo, 'lsCoef') && ~isempty(extraInfo.lsCoef)
                output = extraInfo;
            else
                output.lsCoef = polyfit(Lp, Ls, 1);
                output.ldCoef = polyfit(Lp, Ld, 1);
            end
            
            if iscell(output)
                tmp = cell2mat(output);
                ts = reshape([tmp.lsCoef], 2, []);
                td = reshape([tmp.ldCoef], 2, []);
                ls1 = repmat(ts(1, :), sampNum, 1);
                ls2 = repmat(ts(2, :), sampNum, 1);
                ld1 = repmat(td(1, :), sampNum, 1);
                ld2 = repmat(td(2, :), sampNum, 1);
            else
                ls1 = output.lsCoef(1);
                ls2 = output.lsCoef(2);
                ld1 = output.ldCoef(1);
                ld2 = output.ldCoef(2);
            end
    
            
            
            deltaLs = Ls - ( ls1.*Lp + ls2 ); 
            deltaLd = Ld - ( ld1.*Lp + ld2 );
            x = [Lp; deltaLs; deltaLd];
            
        case 'lpsd'
            x = [Lp; Ls; Ld];
            
        case 'reflectivity'
            D = stpGen1DDiff(sampNum, 1, 1);
            x = [D*Lp; D*Ls; D*Ld];
            
        otherwise
            validatestring(mode, {'Lpsd_fit', 'Lpsd', 'reflectivity'});
    end
    
end