function [f, g] = bsLinearTwoNormWithError(x, data)
    
    z = data.A * x - (data.B - data.y);

    f = sum(z.^2, 1);
    g = 2*(data.A' * z) ;
end