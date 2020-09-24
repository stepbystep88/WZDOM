function res = bsPhase90Shift(data)
% 90 degree phase shift 
    nTrace = size(data, 2);
    
    res = zeros(size(data));
    
    parfor i = 1 : nTrace
        y = hilbert(-data(:, i));
        res(:, i) = imag(y);
    end
    
end