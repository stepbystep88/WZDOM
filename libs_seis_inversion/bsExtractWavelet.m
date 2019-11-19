function [wavelet, GPostInvParam] = bsExtractWavelet(GPostInvParam, timeLine, wellLogs, type)

    
    if ~exist('type', 'var')
        type = 'ricker';
    end
    
    wells = cell2mat(wellLogs);
    wellNum = length(wells);      % the number of wells
    
    inIds = [wells.inline];
    crossIds = [wells.crossline];
    
    sampNum = GPostInvParam.upNum + GPostInvParam.downNum; 
    % horizon of the traces at well location
    horizon = bsCalcHorizonTime(timeLine{GPostInvParam.usedTimeLineId}, inIds, crossIds);
    
    pos = bsCalcT0Pos(GPostInvParam, GPostInvParam.postSeisData.segyInfo, horizon);
    
    % read seismic data
    [postSeisData, GPostInvParam.postSeisData.segyInfo] = bsReadTracesByIds(...
        GPostInvParam.postSeisData.segyFileName, ...
        GPostInvParam.postSeisData.segyInfo, ...
        inIds, ...
        crossIds, ...
        pos, ...
        sampNum-1);
    
    freq = bsGetMainFreq(postSeisData, GPostInvParam.dt);
    GPostInvParam.waveletFreq = freq;
    GPostInvParam.isNormal = 0;
    
    % create wavelet
    switch type
        case 'ricker'
            wave = s_create_wavelet({'type','ricker'}, {'frequencies', freq}, {'step', GPostInvParam.dt}, {'wlength', 80});
        case 'zero-phase'
            wave = s_create_wavelet({'type','zero-phase'}, {'frequencies', freq-20,freq-10,freq+10,freq+10}, {'step', GPostInvParam.dt}, {'wlength', 120});
    end
    wavelet = wave.traces;          
    GPostInvParam.wavelet = wavelet;
    
    % get synthetic data
    scaleFactors = zeros(wellNum, 1);
    similarities = zeros(wellNum, 1);
    
    for i = 1 : wellNum
        synData = bsGetSynthetic(GPostInvParam, wellLogs{i}, horizon(i));
        
%         figure;
%         plot(1:sampNum-1, synData/max(abs(synData)), 'r', 'linewidth', 2); hold on;
%         plot(1:sampNum-1, postSeisData(:, i)/max(abs(postSeisData(:, i))), 'k', 'linewidth', 2);
%         legend('Synthetic', 'Real seismic data');
        
%         scaleFactors(i) = bsCalcScaleFactor(postSeisData(:, i), synData);
        scaleFactors(i) = bsComputeGain(postSeisData(:, i), synData);
        correlation = corrcoef(synData, postSeisData(:, i));
        similarities(i) = correlation(1, 2);
    end
    
    [~, index] = bsMaxK(similarities, ceil(0.8*wellNum));
    meanScaleFactor = mean(scaleFactors(index));

    wavelet = wavelet * meanScaleFactor;
    
%     for i = 1 : wellNum
%         synData = bsGetSynthetic(GPostInvParam, wellLogs{i}, horizon(i));
%         
%         figure;
%         plot(1:sampNum-1, synData*meanScaleFactor, 'r', 'linewidth', 2); hold on;
%         plot(1:sampNum-1, postSeisData(:, i), 'k', 'linewidth', 2);
%         legend('Synthetic', 'Real seismic data');
%         
%     end
    
end


function [B, index] = bsMaxK(A, k, dim)
    
    if ~exist('dim', 'var')
        dim = -1;
    end
        
    [B, I] = sortrows(A, dim);
    B = B(1:k, :);
    index = I(1:k);
end

function freq = bsGetMainFreq(postSeisData, dt)
    N = 512; 
    j = 0 : N - 1; 

    wellNum = size(postSeisData, 2);
    bestft = zeros(wellNum, 1);
    for i = 1 : wellNum
        Fs = 1000.0 / dt; 
        fr = fft(postSeisData(:, i), N);
        mag = sqrt(real(fr).^2 + imag(fr).^2);
        f = j * Fs / N;  
        bestft(i, 1) = min( f(  mag(:, 1) == max(mag)  ) );
    end

    freq = mean( bestft );
    fprintf('The main frequency is :%.2f\n', freq);
        
end
    
function synseis = bsGetSynthetic(GPostInvParam, wellInfo, horizon)
    sampNum = GPostInvParam.upNum + GPostInvParam.downNum; 
    welllog = wellInfo.wellLog;
    dist = horizon - welllog(:, GPostInvParam.indexOfTimeInWellData);
    [~, index] = min(abs(dist));
    s = index - GPostInvParam.upNum;
    trueLog = welllog(s : s+sampNum-1, 1);
    
    model = bsPostPrepareModel(GPostInvParam, wellInfo.inline, wellInfo.crossline, horizon, trueLog, []);
    
    synseis = model.G * model.trueX;
end

function scaleFactor = bsCalcScaleFactor(v1, v2)
    scaleFactor = (v2' * v1) / (v2' * v2);
end

function [waveletGain] = bsComputeGain(post, syn)
    [row, column] = size(syn);
    if column < row
        column = row;
    end
    synDiff=zeros(column-1,1);
    postDiff = zeros(column-1,1);
    for i=2:1:column
        synDiff(i) = syn(i) - syn(i-1);
        postDiff(i) = post(i)-post(i-1);
    end

    synPeaTroValue=zeros(column-1,1);
    postPeaTroValue=zeros(column-1,1);
    k=1;m=1;
    for i=1:1:column-1
        if (synDiff(i) > 0) && (synDiff(i+1) < 0)
            synPeaTroValue(k) = abs(syn(i+1));
            k=k+1;
        end
        if (synDiff(i) < 0) && (synDiff(i+1) > 0)
            synPeaTroValue(k) = abs(syn(i+1));
            k=k+1;
        end
        if (postDiff(i) > 0) &&( postDiff(i+1) < 0)
            postPeaTroValue(k) = abs(post(i+1));
            m=m+1;
        end
        if (postDiff(i) < 0) && (postDiff(i+1) > 0)
            postPeaTroValue(k) = abs(post(i+1));
            m=m+1;
        end
    end

    synPeaTroValue = sort(synPeaTroValue,'descend');
    postPeaTroValue = sort(postPeaTroValue,'descend');

    waveletGain  = mean(postPeaTroValue(postPeaTroValue >0)) / mean(synPeaTroValue(synPeaTroValue >0));
end


    
