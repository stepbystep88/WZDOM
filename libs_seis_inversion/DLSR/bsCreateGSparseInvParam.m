function options = bsCreateGSparseInvParam(DIC, GTrainDICParam, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a struct to control the DLSR-based inversion process
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Dec 2019
% -------------------------------------------------------------------------
    p = inputParser;
    
    if iscell(DIC)
        GTrainDICParam.nAtom = size(DIC{1}, 2);
    else
        GTrainDICParam.nAtom = size(DIC, 2);
    end
    
    
    % the sparsity controlling the sparse representation 
    addParameter(p, 'sparsity', 1); 
    
    % the step size of sliding window
    addParameter(p, 'stride', 1); 
    
    % the way of reconstruction, could be 'equation' or 'simpleAvg'
    % see bsPreInv1DTraceByDLSR.m for details
    addParameter(p, 'reconstructType', 'simpleAvg'); 
    
    % the trained dictionary
    addParameter(p, 'trainDICParam', GTrainDICParam); 
    
    addParameter(p, 'DIC', DIC);
    
	 % the coeficient of normalization, for training CSR dictionary
	addParameter(p, 'rangeCoef', []); 

    addParameter(p, 'output', []);
    
    addParameter(p, 'nMultipleTrace', 1);
    
    addParameter(p, 'nNeibor', 2);
    
    % whether to rebuild the inversion results by using dictionary
    addParameter(p, 'isSparseRebuild', 0); 
    
    % for flag = 'one': train only one dictionary. It mostly used
    % in the poststack inversion
    
    % for flag='ssr': train different dictionaries for different elastic parameters
    
    % for flag='csr': train a joint dictionary which contains the coherence of
    % different elastic parameters
%     validatestring(lower(flag), {'ssr', 'csr', 'one'});
    
    if strcmpi(GTrainDICParam.flag, 'csr')

        % whether to modify the joint dictionary and sparse representation
        % see Eq. 9 Sparse Representation for Color Image Restoration for details
        addParameter(p, 'isModifiedDIC', 0); 
    end
    
    p.parse(varargin{:});  
    options = p.Results;
    options.flag = GTrainDICParam.flag;
    
end

