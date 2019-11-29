function [ref] = bsReflectionProfile(Ips)      

    [sampNum, traceNum] = size(Ips);
    ref = zeros(sampNum-1, traceNum);

    for j = 1 : traceNum
        ref(:,j) = bsReflectionOneTrace(Ips(:,j));
    end
    
    function [ref] = bsReflectionOneTrace(Ip)      

        ref = zeros(sampNum-1, 1);

        for k = 1 : sampNum-1
            ref(k, 1)= (Ip(k+1) - Ip(k)) ./ (Ip(k+1) + Ip(k));
        end

    end
end

