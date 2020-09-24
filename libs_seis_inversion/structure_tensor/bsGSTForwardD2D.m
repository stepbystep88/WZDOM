function [s] = bsGSTForwardD2D(D, r, alpha)
    [n1, n2] = size(r);
    s = zeros(n1, n2);
    
    for i = 1 : n1
        for j = 1 : n2
            e11 = D(i, j, 1) * alpha;
            e12 = D(i, j, 2) * alpha;
            e22 = D(i, j, 3) * alpha;

            r00 = r(i, j);
            if j == 1 && j == 1
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

            ra = r00 - r11;
            rb = r01 - r10;
            r1 = ra - rb;
            r2 = ra + rb;
            s1 = e11 * r1 + e12 * r2;
            s2 = e12 * r1 + e22 * r2;
            sa = s1 + s2;
            sb = s1 - s2;
            s(i, j) = s(i, j) + sa;
            if j > 1
                s(i, j-1) = s(i, j-1) - sb;
            end

            if i > 1
                s(i - 1, j) = s(i - 1, j) + sb;
            end

            if i > 1 && j > 1
                s(i - 1, j - 1) = s(i - 1, j - 1) - sa;
            end

            if i == n1 && j < n2
                s(i, j) = s(i, j) + (sb - sa);
            elseif i < n1 && j == n2
                s(i, j) = s(i, j) - (sa + sb);
            end

        end
    end
end