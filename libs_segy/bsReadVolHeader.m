function [GSegyInfo] = bsReadVolHeader(fileName, GSegyInfo)
%% write volume header into a segy
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
% Input
% fileName          file name of a segy file
% GSegyInfo         a struct containing some basical information of the
% segy file, it will also be output 
%
% Output
% volHeader         volume header
% GSegyInfo         a struct containing some basical information of the segy file
% -------------------------------------------------------------------------
    
    fid = fopen(fileName, 'r', 'ieee-be');
    if fid < 0
        error('The Segy file of name %s opens failed!', fileName);
    end
    
    GSegyInfo.fid = fid;
    GSegyInfo.mode = 'r';
    GSegyInfo.fileName = fileName;
    
    fileInfo = dir(fullfile(fileName));                 % get the basical information of input segy file
    fileSize = fileInfo.bytes;                          % get the size of input segy file
    
    volHeader.baseInfo = fread(GSegyInfo.fid, 3200, 'uchar');     % ASCALL text information
    volHeader.txtInfo = ebcdic2ascii(volHeader.baseInfo);
    
    volFirst = fread(GSegyInfo.fid, 3, 'int32');                  % the first 3 numbers are int32
    volLast = fread(GSegyInfo.fid, (400-12)/2, 'int16');          % the rest numbers are int16
    volBase = [volFirst; volLast];                      % combine them toghter

    volHeader.lastInfo = volBase;
    volHeader.sampInv = volBase(GSegyInfo.sampInvId);                     % the sampling interval
    volHeader.sampNum = volBase(GSegyInfo.sampNumId);                     % the number of samples at each trace
%     volHeader.sampNum = 180;
    
    volHeader.dataForm = volBase(GSegyInfo.dataFormId);                   % data format, 1=IBM, 5=IEEE
%     volHeader.sampNum = 100;
    volHeader.startTime = volBase(GSegyInfo.startTimeId);
    volHeader.traceNum = (fileSize - 3600) / (240 + volHeader.sampNum*4);     % 
    volHeader.sizeTrace = volHeader.sampNum * 4;        % 4 == size of float type 
    
    GSegyInfo.volHeader = volHeader;
end