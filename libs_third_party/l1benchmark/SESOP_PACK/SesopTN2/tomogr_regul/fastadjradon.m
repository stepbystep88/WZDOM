function P = fastadjradon(radontrans, par)
% We need  the angles to be from  [-45,135]
    res = par.res;
    angles = par.angles;
    k = par.k_small_bins;
    Ps = zeros(res);
    Pb = zeros(res); 
    for angleidx = 1:length(angles)
        angle = angles(angleidx);
		
% 		if angle > 135, angle=angle-180; end  % We put the angle to  [-45,135]
% 		if angle > 135, angle=angle-180; end  % We put the angle to  [-45,135]
% 		if angle < -45, angle=angle+180; end  % We put the angle to  [-45,135]
% 		if angle < -45, angle=angle+180; end  % We put the angle to  [-45,135]
		
        if angle > 45
            angle = angle - 90;
            over45 = true;
        else
            over45 = false;
        end
        COSA = cos(angle * pi / 180);
        SINA = sin(angle * pi / 180);
        bigbins = radontrans(:,angleidx);
        smallbins = zeros(floor((length(bigbins)-1)/(COSA/k))+1,1);
        pos = ((length(bigbins)-1)-(length(smallbins)-1)*COSA/k)/2+1;
        bigidx = 1;
        maxbigidx = length(bigbins)-1;
        for i = 1:length(smallbins)
            if (pos-bigidx) < COSA
                if (bigidx+1-pos) < COSA
                    smallbins(i) = bigbins(bigidx)*(COSA-(pos-bigidx))+...
                                   bigbins(bigidx+1)*(COSA-((bigidx+1-pos)));
                else
                    smallbins(i) = bigbins(bigidx)*(COSA-(pos-bigidx));
                end
            else
                smallbins(i) = bigbins(bigidx+1)*(COSA-((bigidx+1-pos)));
            end
            pos = pos + COSA/k;
            bigidx = min(floor(pos),maxbigidx);
        end
        
        smallbins = smallbins / COSA^2; % Normalization

        % go over rows
        os = -((mean([1 res])-1)*COSA-(mean([1 res])-(1:res))*SINA)/(COSA/k); % offset from center
        firstsmallidx = round(mean([1 length(smallbins)]) + os) - 1;
        if over45 == false
            fastbploop2addition(smallbins, Ps, uint32(firstsmallidx), uint32(k));    
        else
            fastbploop2addition(smallbins, Pb, uint32(firstsmallidx), uint32(k));    
        end
    end
    P = Ps' + flipud(Pb);
end