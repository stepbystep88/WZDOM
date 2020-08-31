function denormal_patches = bsSparseDenormalization(normalizationMode, new_patches, input, min_values, max_values)
    switch normalizationMode
        case 'self_mean_sigma'
            denormal_patches = new_patches .* input.sigma_value + input.mean_value;
        case 'whole_data_max_min'
            denormal_patches = new_patches .* (max_values - min_values) + min_values;
        case 'whole_data_mean_sigma'
            denormal_patches = new_patches .* max_values + min_values;  
        case 'feat_mean_sigma'
            denormal_patches = new_patches .* max_values + min_values;
        
        case 'feat_max_min'
%             denormal_patches = new_patches .* (max_value - min_value) + min_value;
            denormal_patches = new_patches .* (max_values - min_values) + min_values;
        otherwise
            denormal_patches = new_patches;
    end
end