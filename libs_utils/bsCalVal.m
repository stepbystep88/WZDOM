function [val_est] = bsCalVal(time_est, time_know, val_know)
    depth = length(time_know);
    if time_est < time_know(1)
        val_est = val_know(1);
    elseif time_est > time_know(depth)
        val_est = val_know(depth);
    else
        d = time_know - time_est;
        n = find(d==0);
        if isempty(n)
            n1 = max(find(d<0));
%             n2 = min(find(d>0));
            n2 = n1 + 1;
            t(1) = time_know(n1);
            t(2) = time_know(n2);
            v(1) = val_know(n1);
            v(2) = val_know(n2);
            val_est = v(1) + (v(2) - v(1))*(time_est - t(1)) / (t(2) - t(1));
        else
            val_est=val_know(n);
        end
    end
end

