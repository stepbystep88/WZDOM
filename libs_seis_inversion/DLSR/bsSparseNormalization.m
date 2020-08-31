function [normal_patches, output] = bsSparseNormalization(normalizationMode, all_patches, min_values, max_values)
    output = [];
    
    switch normalizationMode
        case 'self_mean_sigma'
            mean_value = mean(all_patches, 1);
            sigma_value = var(all_patches, 0, 1);
            
            output.mean_value = repmat(mean_value, size(all_patches, 1), 1);
            output.sigma_value = repmat(sigma_value, size(all_patches, 1), 1);

            normal_patches = (all_patches - output.mean_value) ./ output.sigma_value;
            
        case 'whole_data_max_min'
            normal_patches = (all_patches - min_values) ./ (max_values - min_values);
        case 'whole_data_mean_sigma'
            normal_patches = (all_patches - min_values) ./ max_values;
            
        case 'feat_mean_sigma'

            normal_patches = (all_patches - min_values) ./ max_values;
        
        case 'feat_max_min'
            normal_patches = (all_patches - min_values) ./ (max_values - min_values);
            
        otherwise
            normal_patches = all_patches;
    end
end