function [G] = bsPreBuildGMatrix(mode, vp, vs, angleData, wavelet, extraInfo)
%% Build matrix G for prestack inversion
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------

    sampNum = size(vp, 1);
    angleTrNum = length(angleData);
    
    [c1, c2, c3] = bsAkiSyn(angleData, vp, vs); 
    G = bsGenGByMode(mode, c1, c2, c3, sampNum, angleTrNum, wavelet, extraInfo);
    
end

function [c1, c2, c3] = bsAkiSyn(angleData, vp, vs)

    angleTrNum = length(angleData);
    sampNum = size(vp, 1);
    sampNum = sampNum - 1;
    
    % gamma
    gamma = zeros(sampNum, 1);
    for i = 1 : sampNum
        gamma(i) = ( vs(i+1) + vs(i) ) / ( vp(i+1) + vp(i) );
    end
    
    % meanTheta
    meanTheta = size(sampNum, angleTrNum);
    Theta2 = size(sampNum, angleTrNum);
    for i = 1 : sampNum
        for j = 1 : angleTrNum
            [Theta2(i,j), ~, ~] = snell(angleData(j), vs(i), vp(i), vs(i+1), vp(i+1));
            meanTheta(i,j) = (Theta2(i,j) + angleData(j)) / 2;
        end
    end
    
    
    % c1, c2, c3
    c2 = zeros(sampNum, angleTrNum);
    c3 = zeros(sampNum, angleTrNum);
    
    c1 = 1+(tan(meanTheta)).^2;
    
    for i = 1 : angleTrNum
        c2(:,i) = -8*gamma.^2.*(sin(meanTheta(:,i))).^2;
    end
    for i = 1 : angleTrNum
        c3(:,i) = -(tan(meanTheta(:,i))).^2 + 4*gamma.^2.*(sin(meanTheta(:,i))).^2;
    end
end

function G = bsGenGByMode(mode, c1, c2, c3, sampNum, angleTrNum, wavelet, extraInfo)
    
    W = bsWaveletMatrix(sampNum-1, wavelet);
    cellG = cell(angleTrNum, 3);

    switch lower(mode)
        case 'lpsd_fit'
            D = bsGen1DDiffOperator(sampNum, 1, 1);
            k = extraInfo.lsCoef(1);
            m = extraInfo.ldCoef(1);
            new_c1 = 0.5*c1 + 0.5*k*c2 + 0.5*m*c3;
            new_c2 = 0.5*c2;
            for i = 1 : angleTrNum
                newC1 = diag( new_c1(1:sampNum-1, i) );
                newC2 = diag( new_c2(1:sampNum-1, i) );
                newC3 = diag( c3(1:sampNum-1, i) );

                cellG{i, 1} = W * newC1 * D;
                cellG{i, 2} = W * newC2 * D;
                cellG{i, 3} = W * newC3 * D;
            end
        case 'lpsd'
            D = bsGen1DDiffOperator(sampNum, 1, 1);
            for i = 1 : angleTrNum
                newC1 = diag( 0.5*c1(1:sampNum-1, i) );
                newC2 = diag( 0.5*c2(1:sampNum-1, i) );
                newC3 = diag( 0.5*c3(1:sampNum-1, i) );

                cellG{i, 1} = W * newC1 * D;
                cellG{i, 2} = W * newC2 * D;
                cellG{i, 3} = W * newC3 * D;


            end

        case 'reflectivity'
            for i = 1 : angleTrNum
                newC1 = diag( 0.5*c1(1:sampNum-1, i) );
                newC2 = diag( 0.5*c2(1:sampNum-1, i) );
                newC3 = diag( 0.5*c3(1:sampNum-1, i) );

                cellG{i, 1} = W * newC1;
                cellG{i, 2} = W * newC2;
                cellG{i, 3} = W * newC3;
            end

        otherwise
            
    end
        
    G = sparse( cell2mat(cellG) );
end

