function GSegyInfo = bsWriteVolHeader(fileName, GSegyInfo)
%% write volume header into a segy file
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
% Input
% volHeader         volume header
% fileName          name of the segy file to save the data
% -------------------------------------------------------------------------

    fid = fopen(fileName, 'w', 'ieee-be');
    if fid < 0
        error('The Segy file of name %s opens failed!', fileName);
    end
    
    GSegyInfo.fid = fid;
    GSegyInfo.mode = 'w';
    GSegyInfo.fileName = fileName;
    
    volHeader = GSegyInfo.volHeader;

    % text information of 3200 bytes
    if(fwrite(GSegyInfo.fid, volHeader.baseInfo, 'uchar') ~= 3200)
        disp(ferror(GSegyInfo.fid))
        error('writing failed');
    end

    volHeader.dataForm = 5;
    volHeader.lastInfo(GSegyInfo.dataFormId) = volHeader.dataForm;      % data format:1/IBM;5/IEEE
    volHeader.lastInfo(GSegyInfo.sampInvId) = volHeader.sampInv;        % the sampling interval(3217~3218 in us)
    volHeader.lastInfo(GSegyInfo.sampNumId) = volHeader.sampNum;        % the number of samples(3221~3222)
    volHeader.lastInfo(GSegyInfo.startTimeId) = volHeader.startTime;    % 采样数量(3221~3222字节 样点数/道数=道长)
    
    GSegyInfo.volHeader = volHeader;
    
    fwrite(GSegyInfo.fid, volHeader.lastInfo(1:3), 'ulong');                     % 400字节二进制数区域
    fwrite(GSegyInfo.fid, volHeader.lastInfo(4:end), 'ushort');         
    
end
