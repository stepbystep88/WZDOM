function [output] = bsTrain1DSparseDualDIC(datas, GTrainDICParam)
%% To train a joint dictionary by using K-SVD algorithm
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: March 2020
% -------------------------------------------------------------------------
    
    nWell = length(datas);
    
    num = size(datas{1}, 1);
	index = 1 : GTrainDICParam.stride : num - GTrainDICParam.sizeAtom + 1;
        
    Ph = zeros(GTrainDICParam.sizeAtom, length(index)*nWell);
    PL = zeros(GTrainDICParam.sizeAtom*3, length(index)*nWell);
    counter = 1;
    
    lowData = [];
    for j = 1 : nWell
        lowData = [lowData, datas{j}(:, 1)];
    end
    max_val = max(lowData(:));
    min_val = min(lowData(:));
    
    for j = 1 : nWell
        data = datas{j};
        

        low = data(:, 1);
        low = (low - min_val) / (max_val - min_val);
        high = data(:, 2);
        
        low1 = conv(low, [1,0,-1], 'same'); % the filter is centered and scaled well for s=3
        low2 = conv(low, [1,0,-2,0,1]/2, 'same');

        Eh = high - low;
        
        for k = index
            ee = k+GTrainDICParam.sizeAtom-1;
            
            Ph(:, counter) = Eh(k : ee);
           	PL(:, counter) = [low(k:ee); low1(k:ee); low2(k:ee)];
            
            counter = counter + 1;
        end
    end

    R = PL * PL'; 
    [B, SS] = eig(R); 
    Permute = fliplr(eye(size(R,1))); 
    SS = Permute * SS * Permute; % so that the eigenvalues are sorted descending
    B = B * Permute; 
    energy = cumsum(diag(SS))/sum(diag(SS)); 
    % figure(1); clf; plot(energy)
    pos=find(energy>0.999);
    B = B(:, 1:pos);

    trainData = B' * PL;
    
    params.data = trainData;
    params.Tdata = GTrainDICParam.sparsity;
    params.dictsize = GTrainDICParam.nAtom;
    params.iternum = GTrainDICParam.iterNum;
    params.memusage = 'high';
    
    % using the third-party toolbox to train the dictionary
    [AL, Q] = ksvd(params, GTrainDICParam.isShowIterInfo);
    
    
    
    AH = Ph * Q' * inv(Q * Q');
    output.AL = AL;
    output.AH = AH;
    output.B = B;
    output.max_val = max_val;
    output.min_val = min_val;
end
