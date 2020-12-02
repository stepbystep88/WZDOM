function pN = bsAddNoise(pF, flag, SNR, G, xTrue, dt, options)

    [sampNum, nTrace] = size(pF);
    pN = zeros(sampNum, nTrace);
    
    if ~isfield(options, 'isMeanZero') || options.isMeanZero
        miu = zeros(1, nTrace);
    else
        miu = sign(randn(1, nTrace))*0.4*rand()*max(abs(pF));
    end
    
    sigma = 1;
    
    for i = 1 : nTrace
        
        if ~isempty(G)
            ndata = bsAdd1d(flag, pF(:, i), SNR, miu(i), sigma, xTrue(:, i), G);
        else
            ndata = bsAdd1d(flag, pF(:, i), SNR, miu(i), sigma, [], []);
        end
        
        
        if isfield(options, 'isBandPass') && options.isBandPass
            ndata = bsButtBandPassFilter(ndata, dt, options.lowFreq, options.highFreq);
        end
        
        pN(:, i) = ndata;
    end
    
    if isfield(options, 'isUseHanming') && options.isUseHanming
        NOISE = pN - pF;
        L = 7;
        op = hamming(L);
        ops = sum(sum(op));
        op = op/ops;
    
        NOISE = conv2(NOISE, op, 'same');
    
        pN = pF + NOISE;
    end
    
end

function pN = bsAdd1d(flag, pF, SNR, miu, sigma, xTrue, G)
    switch flag
        case 1
            % non-white white gaussian noise
%             pN = awgn(pF, SNR, 'measured');
            fcn = @(m, n) (normrnd(miu, sigma, m, n));
            pN = bsGenNoiseByFcn(pF, SNR, fcn);
            
        case 2
            % non-white laplacian noise
            fcn = @(m, n) laprnd(m, n, miu, sigma);
            pN = bsGenNoiseByFcn(pF, SNR, fcn);
             
        case 3
            % non-white white gaussian noise
            % gev type 1
            K = rand() * 0.4 + 0.3;
%             K = 0.5;
            sigma = 1;
%             miu = 0;
            fcn = @(m, n) gevrnd(-K, sigma, miu, m, n);
            pN = bsGenNoiseByFcn(pF, SNR, fcn);
        case 4
            sigma = 1;
%             miu = 0;
            fcn = @(m, n) gevrnd(0, sigma, miu, m, n);
            pN = bsGenNoiseByFcn(pF, SNR, fcn);
            
        case 5
            % gev type 1
            K = rand() * 0.4 + 0.3;
%             K = 0.5;
            sigma = 1;
%             miu = 0;
            fcn = @(m, n) gevrnd(K, sigma, miu, m, n);
            pN = bsGenNoiseByFcn(pF, SNR, fcn);
            
        case 6
            % model noise
            xModel = awgn(xTrue, SNR, 'measured');
            xModel = xTrue + 0.1*(xModel - xTrue);
            pN = G * xModel;  
            
            NOISE = pN - pF;
            fcn = @(m, n) (NOISE - mean(NOISE) + miu);
            
            pN = bsGenNoiseByFcn(pF, SNR, fcn);
    end
end

function ndata = bsGenNoiseByFcn(data, SNR, fcn)
    [m, n] = size(data);
    
    NOISE = fcn(m, n);
    NOISE = NOISE - mean(NOISE);
% %     NOISE = NOISE;
    SNR = 10^(SNR / 10);
    
    signal_power = 1/length(data) * sum(data .^ 2);
    noise_power = 1/length(NOISE) * sum(NOISE .^2);
    noise_variance = signal_power / SNR;
    NOISE = sqrt(noise_variance) / sqrt(noise_power) * NOISE;
    ndata = data + NOISE;
    
end