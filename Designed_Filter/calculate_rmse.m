function rmse_val = calculate_rmse(original, processed)
    % calculate_rmse - Calculates the Root Mean Squared Error (RMSE) 
    % between the original and processed images.
    % Inputs:
    %   processed - Processed (noisy or filtered) image.
    %   original  - Original (reference) image.
    % Output:
    %   rmse_val  - Root Mean Squared Error (RMSE).

    mse = mean((original(:) - processed(:)).^2);  % Mean Squared Error
    rmse_val = sqrt(mse);  % Root Mean Squared Error
end

