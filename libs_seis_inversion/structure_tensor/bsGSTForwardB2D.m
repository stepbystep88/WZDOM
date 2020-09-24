function [b] = bsGSTForwardB2D(r)
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