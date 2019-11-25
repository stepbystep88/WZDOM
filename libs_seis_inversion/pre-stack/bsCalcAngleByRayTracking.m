function [meanTheta, angleData] = bsCalcAngleByRayTracking(GPreInvParam, tDepth, tVp, tVs, tRho)
%% ray tracking
    
    offsetInteval = (GPreInvParam.offsetMax - GPreInvParam.offsetMin) / GPreInvParam.newSuperTrNum;
    xoff = zeros(1, GPreInvParam.newSuperTrNum);
    xoff(1) = GPreInvParam.offsetMin + offsetInteval/2;
    for i = 2 : GPreInvParam.newSuperTrNum
        xoff(i) = xoff(i - 1) + offsetInteval;
    end
    sampNum = length(tDepth);                    
    
    if length(tVp) > length(tDepth)
        % assume vp starts from the surface
        num = length(tVp);
        depth = zeros(num, 1);
        for t = 1 : num - 1
            depth(t+1) = depth(t) + tVp(t) * GPreInvParam.dt * 0.001 * 0.5;
        end
        [~, index] = min( abs(depth - tDepth(1) ) );
        
        z0Num = index - 1;
        vp = tVp(1:z0Num+sampNum, 1);
        depth = depth(1:z0Num+sampNum, 1);
        vs = [ones(z0Num, 1) * tVs(1); tVs];
        rho = [ones(z0Num, 1) * tRho(1); tRho];
    else
        dz = tVp(1) * 0.001 * GPreInvParam.dt * 0.5;            % 
        z0 = 0 : dz : tDepth(1);                                % 
        z0Num = length(z0);                                     % 

        vp0 = zeros(1, z0Num); vs0 = zeros(1, z0Num); rho0 = zeros(1, z0Num);
        vp0(:) = tVp(1); vs0(:) = tVs(1); rho0(:) = tRho(1);
        depth = [z0'; tDepth]; vp = [vp0'; tVp]; vs = [vs0'; tVs]; rho = [rho0'; tRho];
    end
    

    startZ = z0Num;    
    ts = zeros(sampNum, GPreInvParam.newSuperTrNum);             
    ps = zeros(sampNum, GPreInvParam.newSuperTrNum);             

    zsrc=0; zrec=0; caprad=1e8; itermax=50; pfan=-1; optflag=1; pflag=1; dflag=2; 

    for j = 1 : sampNum
        zd = depth(startZ + j);                      
        [t, p] = traceray_pp(vp, depth, zsrc, zrec, zd, xoff, caprad, pfan, itermax, optflag, pflag, dflag);
        ts(j, :) = t;
        ps(j, :) = p;
    end

     %Zoepritz equation
    [meanTheta] = Zoep_syn(xoff, tDepth, startZ, ps, vp, vs, rho);   
    
    %% angle should be smaller than max angle
    setMaxAngle = GPreInvParam.maxAngle * pi / 180;
    maxAngle = max(max(meanTheta));
    minAngle = min(min(meanTheta));
    if( maxAngle > setMaxAngle)
        maxAngle = setMaxAngle;
    end
    
    angleScale = minAngle : ( maxAngle - minAngle)/(GPreInvParam.angleTrNum) : maxAngle;
    angleData = zeros(1, GPreInvParam.angleTrNum);
    for i = 1: 1 : (length(angleScale)-1)
        angleData(i) = (angleScale(i) + angleScale(i+1)) / 2;
    end

end