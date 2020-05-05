function [data, B] = bsReduction(data, pos)

    R = data * data'; 
    [B, SS] = eig(R); 
    Permute = fliplr(eye(size(R, 1))); 
    SS = Permute * SS * Permute; % so that the eigenvalues are sorted descending
    B = B * Permute; 
    energy = cumsum(diag(SS))/sum(diag(SS)); 
    % figure(1); clf; plot(energy)
    
    if nargin == 2 && isempty(pos)
        pos=find(energy>0.999, 1);
    end
%     pos = 15;
    B = B(:, 1:pos);
    data = B' * data;
        
end