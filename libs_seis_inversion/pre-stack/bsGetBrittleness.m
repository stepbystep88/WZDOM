function data = bsGetBrittleness(vp, vs, rho)
    wIp = vp .* rho;
    wIs = vs .* rho;
    wrho = rho;
    wPossion = (wIp.^2 - 2*wIs.^2) ./ (2*wIp.^2 - wIs.^2);
    wPossion(wPossion<0, 1) = 0;
    
    minPossin = min(wPossion(:));
    maxPossin = max(wPossion(:));
    
    %% 杨氏模量计算公式  
    wYangShi = wIs.^2 .* (3*wIp.^2 - 4*wIs.^2) ./ (wrho.*(wIp.^2 - wIs.^2));
    minYangShi = min(wYangShi(:));
    maxYangShi = max(wYangShi(:));
    
    %% 计算脆性
    wCuiXing = 50 * ( (wYangShi - minYangShi)/(maxYangShi - minYangShi) + (wPossion-maxPossin)/(minPossin-maxPossin));

    data = wCuiXing;
    
end