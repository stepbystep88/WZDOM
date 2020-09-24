function xk = bsGSTSmooth2D(D, f, alpha, maxIter)
%% 利用结构张量光滑数据
    if nargin < 3
        alpha = 10;
    end
    
    if nargin < 4
        maxIter = 30;
    end
    
    % 求b
%     b = BTBf(f);
%     b = medfilt2(f, [4 4]);
    b = f;
    
%     xk = bsCGMultiTraces(f, b, @forwardAx, maxIter);
    xk = b;
    rk = b - forwardAx(D, xk, alpha);
    pk = rk;
    
    for k = 1 : maxIter
        Apk = forwardAx(D, pk, alpha);
        alphak = sum(rk .^ 2) / sum(pk .* Apk);
%         alphak = 0.01;
        xk = xk + alphak * pk;
        rk_new = rk - alphak * Apk;
        
        betak = sum(rk_new.^2) / sum(rk.^2);
        pk = rk_new + betak * pk;
        rk = rk_new;
    end
    
end

function b = BTBf(r)
    [n1, n2] = size(r);
    b = zeros(n1, n2);
    
    for i = 1 : n1
        for j = 1 : n2
            
            r00 = r(i, j);
            if j == 1
                r01 = r(i, 1);
                
                if i == 1
                    r10 = r(1, j);
                    r11 = r(1, 1);
                else
                    r10 = r(i - 1, j);
                    r11 = r(i - 1, 1);
                end
            
            else
                r01 = r(i, j - 1);
                
                if i == 1
                    r10 = r(1, j);
                    r11 = r(1, j - 1);
                else
                    r10 = r(i - 1, j);
                    r11 = r(i - 1, j - 1);
                end
            end
            
            b(i, j) = (r00 + r01 + r10 + r11);
        end
    end
    
end

function s = forwardAx(D, r, alpha)
    
    s = bsGSTForwardD2D(D, r, alpha) + r;
end