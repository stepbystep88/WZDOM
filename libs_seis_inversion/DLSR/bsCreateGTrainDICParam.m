function options = bsCreateGTrainDICParam(flag, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a struct to control the dictionary training process
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Dec 2019
% -------------------------------------------------------------------------
    if nargin < 1
        flag = 'ssr';
    end
    
    p = inputParser;
    
    % the length of a atom (column vector of the dictionary matrix)
    addParameter(p, 'sizeAtom', 50); 
    
    % the number of atoms of the dictionary
    addParameter(p, 'nAtom', 2000); 
    
    % the sparsity controlling the sparse representation in the training
    % process
    addParameter(p, 'sparsity', 3); 
    
    % the step size of sliding window
    addParameter(p, 'stride', 1); 
    
    % the number of iterations to train a dictionary
    addParameter(p, 'iterNum', 10); 
    
    % the coeficient of a low-pass filter which is used to pre-process the
    % well-log data
    addParameter(p, 'filtCoef', 0.4); 
    
    % whether to show the iteration information
    addParameter(p, 'isShowIterInfo', 0); 
    
    addParameter(p, 'isRebuild', 0); 
        
    % indicate the name of the trained dictionary
    addParameter(p, 'title', []);
    
    addParameter(p, 'isAddLocInfo', 0);
    
    addParameter(p, 'isAddTimeInfo', 0);
    
    switch lower(flag)
        case 'one'
            % for flag = 'one': train only one dictionary. It mostly used
            % in the poststack inversion
            addParameter(p, 'dicSavePath', './TrainedDictionaries/poststack/'); 
            addParameter(p, 'normalizationMode', 'none'); 
        case 'ssr'
            % for flag='ssr': train different dictionaries for different elastic parameters
            % the path of saving the trained dictionary
            addParameter(p, 'dicSavePath', './TrainedDictionaries/ssr/'); 
            addParameter(p, 'normalizationMode', 'none'); 
            
        case 'csr'
            % for flag='csr': train a joint dictionary which contains the coherence of
            % different elastic parameters
            addParameter(p, 'dicSavePath', './TrainedDictionaries/csr/'); 
            addParameter(p, 'normalizationMode', 'whole_data_max_min'); 
            addParameter(p, 'feature_reduction', 0); 
            
        case 'dual'
            addParameter(p, 'dicSavePath', './TrainedDictionaries/dual/'); 
            addParameter(p, 'normalizationMode', 'none'); 
        otherwise
            validatestring(flag, {'ssr', 'csr', 'one', 'dual'});
    end
    
    p.parse(varargin{:});  
    options = p.Results;
    options.flag = flag;
    
    
end

