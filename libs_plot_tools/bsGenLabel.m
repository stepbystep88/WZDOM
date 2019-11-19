function [label, data] = bsGenLabel(minVal, maxVal, sampNum, num)
    %% Éú³É¿Ì¶È
%     l = linspace(0, sampNum, sampNum+1);
%     d = linspace(minVal, maxVal, sampNum+1);
    
    step = sampNum / num;
    
    label = zeros(1, num);
    data = zeros(1, num);
    
    index = sampNum+1;
    
    x = 0;
    for i = 1 : num
%         label(num-i+1) = l(index);
%         data(num-i+1) = round(d(index));
%         index = index - step;
        label(i) = round(x);
        
        data(i) = round(round(x) / sampNum * (maxVal - minVal) + minVal);
        x = x + step;
    end
    
end