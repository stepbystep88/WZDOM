function [range, scale, dataIndex, dataColorTbl, attName] = ...
    bsGetInfoByType(GShowProfileParam, GInvParam, type)

    dataIndex = [];
    scale = 1;
    dataColorTbl = [];
        
    switch lower(type)
        case 'ip'

            range = GShowProfileParam.range.ip;
            scale = 1000;
            dataIndex = GInvParam.indexInWellData.ip;
            dataColorTbl = GShowProfileParam.colormap.ip;
            attName = 'Impedance (g/cm^3\cdotkm/s)';

        case 'vp'
            range = GShowProfileParam.range.vp;
            scale = 1000;
            dataIndex = GInvParam.indexInWellData.vp;
            dataColorTbl = GShowProfileParam.colormap.vp;
            attName = 'V_P km/s';

        case 'vs'
            range = GShowProfileParam.range.vs;
            scale = 1000;
            dataIndex = GInvParam.indexInWellData.vs;
            dataColorTbl = GShowProfileParam.colormap.vs;
            attName = 'V_S km/s';

        case {'rho', 'density'}
            range = GShowProfileParam.range.rho;
            dataIndex = GInvParam.indexInWellData.rho;
            dataColorTbl = GShowProfileParam.colormap.rho;
            attName = 'Density g/cm^3';

        case {'vp_vs', 'vpvs_ratio'}
            attName = 'Vp/Vs ratio';
            range = GShowProfileParam.range.vpvs_ratio;
            dataIndex = GInvParam.indexInWellData.vpvs_ratio;
            dataColorTbl = GShowProfileParam.colormap.vpvs_ratio;

        case {'possion'}
            attName = 'Possion ratio';
            range = GShowProfileParam.range.possion;
            dataIndex = GInvParam.indexInWellData.possion;
            dataColorTbl = GShowProfileParam.colormap.possion;
        case 'seismic'
            attName = 'Seismic (Amplitude)';
%                 [wellPos, ~] = bsFindWellLocation(wellLogs, profile.inIds, profile.crossIds);
            range = GShowProfileParam.range.seismic;
            dataColorTbl = GShowProfileParam.colormap.seismic;

        otherwise
            validatestring({'seismic', 'ip', ...
                'vp', 'vs', 'rho', 'density', ...
                'vp_vs', 'vpvs_ratio', 'possion'});
    end
end