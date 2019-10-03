function [z] = bsLevy(n, m, beta)
% This function implements Levy's flight. 
% For more information see 
%'Multiobjective cuckoo search for design optimization Xin-She Yang, Suash Deb'. 
% Coded by Hemanth Manjunatha on Nov 13 2015.
% Input parameters
% m     -> Number of steps 
% n     -> Number of Dimensions 
% beta  -> Power law index  % Note: 1 < beta < 2
% Output 
% z     -> 'n' levy steps in 'm' dimension

    num = gamma(1+beta).*sin(pi*beta/2); % used for Numerator 
    den = gamma((1+beta)/2).*beta.*2.^((beta-1)/2); % used for Denominator
    sigma_u = (num./den) .^ (1./beta);% Standard deviation
    
    if length(beta) ~= 1
        sigma_u = repmat(sigma_u, n, 1);
    end

    u = randn(n, m) .* sigma_u;
    v = randn(n, m);
    
    z = u./(abs(v).^(1./beta));
    
end
