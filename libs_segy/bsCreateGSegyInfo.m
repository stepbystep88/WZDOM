function options = bsCreateGSegyInfo(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a struct to save the information of a segy file 
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Dec 2019
% -------------------------------------------------------------------------
    p = inputParser;
    
    % the index (in line header) of sample interval
    addParameter(p, 'sampInvId', 6); 
    
    % the index (in line header) of number of samples in each trace
    addParameter(p, 'sampNumId', 8);  
    
    % the index (in line header) of data format (IEEE or IBM)
    addParameter(p, 'dataFormId', 10);
    
    % the index (in line header) of start time of the sgy file, this entry
    % has not been used actually
    addParameter(p, 'startTimeId', 50);
    
    % the start byte index (in trace header) of inline 
    addParameter(p, 'inlineId', 9);
    
    % the start byte index (in trace header) of crossline 
    addParameter(p, 'crosslineId', 21);
    
    % the start byte index (in trace header) of x coordinate
    addParameter(p, 'xCoordId', 181);
    
    % the start byte index (in trace header) of y coordinate
    addParameter(p, 'yCoordId', 185);
    
    % the start byte index (in trace header) of offset (for prestack or angle gather)
    addParameter(p, 'offsetId', 37);
    
    % the start byte index of offset
    addParameter(p, 'traceSampNumId', 113);
    
    % whether to save the corrdinate information when write a segy file
    addParameter(p, 'isSaveCoordInfo', 1);
    
    % whether read the data as -data
    addParameter(p, 'isNegative', 0);
    
    % the line header information of a segyfile
    addParameter(p, 'volHeader', []);
    
    % the trace header information of a segyfile
%     addParameter(p, 'traceHeader', []);

    % the file handle of a segy file
    addParameter(p, 'fid', -1);
    
    % the file name of a segy file
    addParameter(p, 'fileName', '');
    
    % the start time of a segy file, must be set by users
    addParameter(p, 't0', 0);
    
    % whether to read the traces from begining, if it is set to 1, the
    % program will ignore the start time information, and just read the
    % data of each trace from the first sample point
    % see bsGetT0Pos.m for details
    addParameter(p, 'isPosZero', 0);

    p.parse(varargin{:});  
    options = p.Results;
    
end

