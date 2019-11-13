function bsShowInvProfiles(GPostInvParam, GPlotParam, GShowProfileParam, profiles)
%% Show the inversion results
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
    
    nProfile = length(profiles);
    [~, traceNum] = size(profiles{1}.data);
    
    figure;
    switch nProfile
        case 1
            set(gcf, 'position', [432   246   981   385]);
            nRow = 1;
            nCol = 1;
            
        case 2
            set(gcf, 'position', [0   0   981   700]);
            nRow = 2;
            nCol = 1;
            
        case 3
            set(gcf, 'position', [0   0   981   700]);
            nRow = 3;
            nCol = 1;
        case 4
            set(gcf, 'position', [0   0   1800   700]);
            nRow = 2;
            nCol = 2;
        case 5
            set(gcf, 'position', [0   0   1800   800]);
            nRow = 3;
            nCol = 2;
        case 6
            set(gcf, 'position', [0   0   1800   800]);
            nRow = 3;
            nCol = 2;
    end
    
        
    for iProfile = 1 : nProfile
        
        profile = profiles{iProfile};
        
        subplot(nRow, nCol, iProfile);
        
        bsShowOneProfile(profile);
        
        
    end
    
    function bsShowOneProfile(profile)
        
        profileData = profile.data;
        if profile.inIds(1) == profile.inIds(2)
            traceIds = profile.crossIds;
        else
            traceIds = profile.inIds;
        end
        
        [sampNum, ~] = size(profileData);
        
        % data preparation
        switch profile.type
            case 'IP'
                
                profileData(profileData<=0) = inf;
                
                if(max(max(profileData)) > 1000)
                    profileData = profileData / 1000;
                end
                
                
                
                range = GShowProfileParam.rangeIP;
                if range(1) > 1000
                    range = range / 1000;
                end
        end
        
        bsShowSubProfile(GPostInvParam, GPlotParam, GShowProfileParam, ...
            profileData, profile.horizon, traceIds(1), ... 
            0, GShowProfileParam.showProfileFiltCoef, profile.name, range, GShowProfileParam.dataColorTbl, [], []);
        
        switch profile.type
            case 'IP'
                
                ylabel(colorbar('fontsize', GPlotParam.fontsize,'fontweight', GPlotParam.fontweight, 'fontname', GPlotParam.fontname), ...
                    'Impedance (g/cm^3\cdotkm/s)', 'fontsize', GPlotParam.fontsize,'fontweight', GPlotParam.fontweight, 'fontname', GPlotParam.fontname);
            case 'Seismic'
                ylabel(colorbar('fontsize', GPlotParam.fontsize,'fontweight', GPlotParam.fontweight, 'fontname', GPlotParam.fontname), ...
                    'Seismic (Amplitude)', 'fontsize', GPlotParam.fontsize,'fontweight', GPlotParam.fontweight, 'fontname', GPlotParam.fontname);
        end
    end
end
