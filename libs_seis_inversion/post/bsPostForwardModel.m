function [GInvParam, poststackFreeNoise, poststackNoise, wavelet, options] ...
    = bsPostForwardModel(GInvParam, trueModel, varargin)

     p = inputParser;
    
    addParameter(p, 'title', '');
    addParameter(p, 'dstPath', sprintf('%s/', GInvParam.modelSavePath));
    addParameter(p, 'noistFlag', 1);
    addParameter(p, 'SNR', 10);
    addParameter(p, 'isRebuild', 0);
    
    addParameter(p, 'waveletType', 'ricker');
    
    addParameter(p, 'isMeanZero', '0');
    addParameter(p, 'isBandPass', false);
    addParameter(p, 'lowFreq', 5);
    addParameter(p, 'highFreq', 100);
    
    addParameter(p, 'isUseHanming', 1);
    
    p.parse(varargin{:});  
    options = p.Results;
    
    fileName = bsGetFileName();
    
    if( exist(fileName, 'file') && ~options.isRebuild)
        matFile = load( fileName, ...
            'wavelet', 'poststackFreeNoise', 'poststackNoise');
        wavelet = matFile.wavelet;
        GInvParam.wavelet = wavelet;
        
        poststackFreeNoise = matFile.poststackFreeNoise;
        poststackNoise = matFile.poststackNoise;
    else
        sampNum = size(trueModel, 1);
        logModel = log(trueModel);
        wavelet = bsGenWavelet(options.waveletType, ...
             GInvParam.waveletFreq, GInvParam.dt, []);
        GInvParam.wavelet = wavelet;
        
        G = bsPostGenGMatrix(GInvParam.wavelet, sampNum);
        poststackFreeNoise = G * logModel;
        
        if( options.SNR == 0)
            poststackNoise = poststackFreeNoise;
        else
            poststackNoise = bsAddNoise(poststackFreeNoise, ...
                options.noistFlag, options.SNR, G, logModel, GInvParam.dt, options);
        end
        
        save(fileName, 'wavelet', 'poststackFreeNoise', 'poststackNoise');
    end
    
    function fileName = bsGetFileName()
        fileName = sprintf('%s/SynModel_MainFreq%.1f_dt_%d_nflag_%d_SNR_%.1f.mat', ...
            options.dstPath, GInvParam.waveletFreq, GInvParam.dt, options.noistFlag, options.SNR);
    end
end