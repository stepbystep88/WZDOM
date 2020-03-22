function bigbins = fastfp(P, angle, k)
    if angle <= 45
        flipbigbins = false;
    elseif angle < 90
        P = P';
        angle = 90 - angle;
        flipbigbins = true;
    elseif angle < 135
        P = rot90(P,-1);
        angle = angle - 90;
        flipbigbins = false;
    else
        P = flipud(P);
        angle = 180 - angle;
        flipbigbins = true;
    end
    COSA = cos(angle * pi / 180);
    SINA = sin(angle * pi / 180);
    res = size(P,1);
    bigbins = zeros(1,ceil(sqrt(2)*(res-1))+1);
    smallbins = zeros(1,floor((length(bigbins)-1)/(COSA/k))+1);
    
    constskips = k*((1:res)-1);
    % go over rows
    os = -((mean([1 res])-1)*COSA-(mean([1 res])-(1:res))*SINA)/(COSA/k); % offset from center
    firstsmallidx = round(mean([1 length(smallbins)]) + os);
    for r = 1:res;
        smallbins(firstsmallidx(r)+constskips) = smallbins(firstsmallidx(r)+constskips) + P(r,:);
    end

    pos = ((length(bigbins)-1)-(length(smallbins)-1)*COSA/k)/2+1;
    bigidx = 1;
    maxbigidx = length(bigbins)-1;
    for i = 1:length(smallbins)
        bigbins(bigidx) = bigbins(bigidx) + smallbins(i)*(1-(pos-bigidx));
        bigbins(bigidx+1) = bigbins(bigidx+1) + smallbins(i)*(pos-bigidx);
        pos = pos + COSA/k;
        bigidx = min(floor(pos),maxbigidx);
    end
    if flipbigbins
        bigbins = fliplr(bigbins);
    end
end