function [GInvParam, prestackFreeNoise, prestackNoise, options] ...
    = bsPreForwardModel(GInvParam, depth, vp, vs, rho, varargin)

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
            'wavelet', 'prestackFreeNoise', 'prestackNoise', 'GInvParam');
        wavelet = matFile.wavelet;
        GInvParam.wavelet = wavelet;
        
        prestackNoise = matFile.prestackNoise;
        prestackFreeNoise = matFile.prestackFreeNoise;
        GInvParam = matFile.GInvParam;
        
    else
        [sampNum, traceNum] = size(vp);
        angleData = linspace(0, GInvParam.maxAngle, GInvParam.angleTrNum + 1);
        angleData = angleData(2:end);
        if angleData(end) > 10
            angleData = angleData / 180 * pi;
        end
            
        GInvParam.angleData = angleData;
        wavelet = bsGenWavelet(options.waveletType, ...
             GInvParam.waveletFreq, GInvParam.dt, []);
        GInvParam.wavelet = wavelet;
        
        prestackFreeNoise = cell(1, traceNum);
        prestackNoise = cell(1, traceNum);
        
        parfor i = 1 : traceNum
            wellData = [depth(:, i), vp(:, i), vs(:, i), rho(:, i)];
        
            % build model parameter
            [x, lsdCoef] = bsPreBuildModelParam(wellData, GInvParam.mode, GInvParam.lsdCoef);

            if GInvParam.isInitDeltaZero
                x(sampNum+1:end) = 0;
            end

            % build forward matrix G
            G = bsPreBuildGMatrix(...
                        GInvParam.mode, ...
                        wellData(:, 2), ...
                        wellData(:, 3), ...
                        angleData, ...
                        GInvParam.wavelet, ...
                        lsdCoef);

            % reshape angle seismic data as a vector
            d = G * x;
            angleSeisData = reshape(d, [], GInvParam.angleTrNum);
            prestackFreeNoise{i} = angleSeisData;
            
            
            
            if( options.SNR == 0)
                prestackNoise{i} = prestackFreeNoise{i};
            else
                prestackNoise{i} = bsAddNoise(prestackFreeNoise{i}, ...
                    options.noistFlag, options.SNR, G, x, GInvParam.dt, options);
            end
        
        end
        
        save(fileName, 'GInvParam', 'wavelet', 'prestackFreeNoise', 'prestackNoise');
        
    end
    
    function fileName = bsGetFileName()
        if ~exist(options.dstPath, 'dir')
            mkdir(options.dstPath);
        end
        
        fileName = sprintf('%s/SynPreModel_MainFreq%.1f_dt_%d_nflag_%d_SNR_%.1f.mat', ...
            options.dstPath, GInvParam.waveletFreq, GInvParam.dt, options.noistFlag, options.SNR);
    end
end