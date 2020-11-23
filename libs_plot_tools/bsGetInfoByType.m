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
            
            if strcmpi(GShowProfileParam.language, 'en')
                attName = 'Impedance (g/cm^3\cdotkm/s)';
            else
                attName = '×è¿¹ (g/cm^3\cdotkm/s)';
            end

        case 'vp'
            range = GShowProfileParam.range.vp;
            scale = 1000;
            dataIndex = GInvParam.indexInWellData.vp;
            dataColorTbl = GShowProfileParam.colormap.vp;
            if strcmpi(GShowProfileParam.language, 'en')
                attName = 'V_P km/s';
            else
                attName = '×Ý²¨ËÙ¶È (km/s)';
            end

        case 'vs'
            range = GShowProfileParam.range.vs;
            scale = 1000;
            dataIndex = GInvParam.indexInWellData.vs;
            dataColorTbl = GShowProfileParam.colormap.vs;
            if strcmpi(GShowProfileParam.language, 'en')
                attName = 'V_S km/s';
            else
                attName = 'ºá²¨ËÙ¶È (km/s)';
            end

        case {'rho', 'density'}
            range = GShowProfileParam.range.rho;
            dataIndex = GInvParam.indexInWellData.rho;
            dataColorTbl = GShowProfileParam.colormap.rho;
            
            if strcmpi(GShowProfileParam.language, 'en')
                attName = 'Density g/cm^3';
            else
                attName = 'ÃÜ¶È (g/cm^3)';
            end

        case {'vp_vs', 'vpvs_ratio'}
            if strcmpi(GShowProfileParam.language, 'en')
                attName = 'Vp/Vs ratio';
            else
                attName = '×Ýºá²¨±ÈÂÊ';
            end
            
            range = GShowProfileParam.range.vpvs_ratio;
            dataIndex = GInvParam.indexInWellData.vpvs_ratio;
            dataColorTbl = GShowProfileParam.colormap.vpvs_ratio;
    
        case {'possion'}
            if strcmpi(GShowProfileParam.language, 'en')
                attName = 'Possion ratio';
            else
                attName = '²´ËÉ±È';
            end
            
            range = GShowProfileParam.range.possion;
            dataIndex = GInvParam.indexInWellData.possion;
            dataColorTbl = GShowProfileParam.colormap.possion;
        case 'seismic'
            if strcmpi(GShowProfileParam.language, 'en')
                attName = 'Seismic (Amplitude)';
            else
                attName = 'µØÕð (Õñ·ù)';
            end
%                 [wellPos, ~] = bsFindWellLocation(wellLogs, profile.inIds, profile.crossIds);
            range = GShowProfileParam.range.seismic;
            dataColorTbl = GShowProfileParam.colormap.seismic;

        otherwise
%             validatestring({'seismic', 'ip', ...
%                 'vp', 'vs', 'rho', 'density', ...
%                 'vp_vs', 'vpvs_ratio', 'possion'});]
            attName = type;
            range = GShowProfileParam.range.other;
            dataColorTbl = GShowProfileParam.colormap.other;
    end
end