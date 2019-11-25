function G = bsPostGenGMatrix(wavelet, sampNum)
%% create forward matrix G for poststack inversion
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
    D = bsGen1DDiffOperator(sampNum, 1, 1);
    W = bsWaveletMatrix(sampNum-1, wavelet);
    G = 0.5 * W * D;
end