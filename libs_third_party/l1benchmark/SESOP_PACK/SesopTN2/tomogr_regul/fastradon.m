function radontrans = fastradon(P, par)
% We need  the angles to be from  [-45,135]
    res = par.res;
    angles = par.angles;
    k = par.k_small_bins;
    nbigbins = ceil(sqrt(2)*(res-1))+1;
    radontrans = zeros(nbigbins,length(angles));
    Ps = P';
    Pb = flipud(P);
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
        smallbins = zeros(floor((nbigbins-1)/(COSA/k))+1,1);
        bigbins = zeros(nbigbins,1);

        % go over rows
        os = -((mean([1 res])-1)*COSA-(mean([1 res])-(1:res))*SINA)/(COSA/k); % offset from center
        firstsmallidx = round(mean([1 length(smallbins)]) + os) - 1;
        if over45 == false
            fastfploop1(smallbins, Ps, uint32(firstsmallidx), uint32(k));
        else
            fastfploop1(smallbins, Pb, uint32(firstsmallidx), uint32(k));
        end

        pos = ((length(bigbins)-1)-(length(smallbins)-1)*COSA/k)/2+1;
        bigidx = 1;
        maxbigidx = length(bigbins)-1;
        for i = 1:length(smallbins)
            if (pos-bigidx) < COSA
                bigbins(bigidx) = bigbins(bigidx) + smallbins(i)*(COSA-(pos-bigidx));
            end
            if (bigidx+1-pos) < COSA
                bigbins(bigidx+1) = bigbins(bigidx+1) + smallbins(i)*(COSA-((bigidx+1-pos)));
            end
            pos = pos + COSA/k;
            bigidx = min(floor(pos),maxbigidx);
        end
        radontrans(:,angleidx) = bigbins / COSA^2; % Normalization
    end
end