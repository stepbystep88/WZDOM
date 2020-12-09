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
                attName = '�迹 \fontname{Times New Roman}(g/cm^3\cdotkm/s)';
            end

        case 'vp'
            range = GShowProfileParam.range.vp;
            scale = 1000;
            dataIndex = GInvParam.indexInWellData.vp;
            dataColorTbl = GShowProfileParam.colormap.vp;
            if strcmpi(GShowProfileParam.language, 'en')
                attName = 'V_P km/s';
            else
                attName = '�ݲ��ٶ� \fontname{Times New Roman}(km/s)';
            end

        case 'vs'
            range = GShowProfileParam.range.vs;
            scale = 1000;
            dataIndex = GInvParam.indexInWellData.vs;
            dataColorTbl = GShowProfileParam.colormap.vs;
            if strcmpi(GShowProfileParam.language, 'en')
                attName = 'V_S km/s';
            else
                attName = '�Შ�ٶ� \fontname{Times New Roman}(km/s)';
            end

        case {'rho', 'density'}
            range = GShowProfileParam.range.rho;
            dataIndex = GInvParam.indexInWellData.rho;
            dataColorTbl = GShowProfileParam.colormap.rho;
            
            if strcmpi(GShowProfileParam.language, 'en')
                attName = 'Density g/cm^3';
            else
                attName = '�ܶ� \fontname{Times New Roman}(g/cm^3)';
            end

        case {'vp_vs', 'vpvs_ratio'}
            if strcmpi(GShowProfileParam.language, 'en')
                attName = 'Vp/Vs ratio';
            else
                attName = '�ݺᲨ����';
            end
            
            range = GShowProfileParam.range.vpvs_ratio;
            dataIndex = GInvParam.indexInWellData.vpvs_ratio;
            dataColorTbl = GShowProfileParam.colormap.vpvs_ratio;
    
        case {'possion'}
            if strcmpi(GShowProfileParam.language, 'en')
                attName = 'Possion ratio';
            else
                attName = '���ɱ�';
            end
            
            range = GShowProfileParam.range.possion;
            dataIndex = GInvParam.indexInWellData.possion;
            dataColorTbl = GShowProfileParam.colormap.possion;
        
        case {'brittleness', 'cuixing', 'cui_xing'}
            if strcmpi(GShowProfileParam.language, 'en')
                attName = 'Brittleness';
            else
                attName = '����';
            end
            
            range = GShowProfileParam.range.brittleness;
            dataIndex = GInvParam.indexInWellData.brittleness;
            dataColorTbl = GShowProfileParam.colormap.brittleness;
        
        case {'toc'}
            if strcmpi(GShowProfileParam.language, 'en')
                attName = 'TOC';
            else
                attName = '�л���̼��';
            end
            
            range = GShowProfileParam.range.toc;
            dataIndex = GInvParam.indexInWellData.toc;
            dataColorTbl = GShowProfileParam.colormap.toc;
            
        case 'seismic'
            if strcmpi(GShowProfileParam.language, 'en')
                attName = 'Seismic (Amplitude)';
            else
                attName = '���� (���)';
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