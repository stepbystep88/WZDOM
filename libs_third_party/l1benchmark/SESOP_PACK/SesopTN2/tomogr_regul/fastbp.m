function P = fastbp(bigbins, angle, k)
    if angle <= 45
        flipP = 0; % Do not flip
    elseif angle < 90
        bigbins = fliplr(bigbins);
        flipP = 1; % Transpose
        angle = 90 - angle;
    elseif angle < 135
        flipP = 2; % Rotate 90 degrees
        angle = angle - 90;
    else
        bigbins = fliplr(bigbins);
        flipP = 3; % Flip up down
        angle = 180 - angle;
    end
    COSA = cos(angle * pi / 180);
    SINA = sin(angle * pi / 180);
    res = floor((length(bigbins)-1)/sqrt(2))+1;
    P = zeros(res);
    smallbins = zeros(1,floor((length(bigbins)-1)/(COSA/k))+1);

    pos = ((length(bigbins)-1)-(length(smallbins)-1)*COSA/k)/2+1;
    bigidx = 1;
    maxbigidx = length(bigbins)-1;
    for i = 1:length(smallbins)
        smallbins(i) = bigbins(bigidx)*(1-(pos-bigidx))+bigbins(bigidx+1)*(pos-bigidx);
        pos = pos + COSA/k;
        bigidx = min(floor(pos),maxbigidx);
    end   
    
    constskips = k*((1:res)-1);
    % go over rows
    os = -((mean([1 res])-1)*COSA-(mean([1 res])-(1:res))*SINA)/(COSA/k); % offset from center
    firstsmallidx = round(mean([1 length(smallbins)]) + os);
    for r = 1:res;
        P(r,:) = smallbins(firstsmallidx(r)+constskips);
    end    
    
    switch flipP
        case 1
            P = P';
        case 2
            P = rot90(P);
        case 3
            P = flipud(P);
    end
end