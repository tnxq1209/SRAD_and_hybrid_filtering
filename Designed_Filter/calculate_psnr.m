% ----- PSNR Calculation Function -----
function psnr_value = calculate_psnr(original, processed)
    % Calculate Peak Signal-to-Noise Ratio (PSNR)
    % Inputs:
    %   original  - Original image (grayscale, normalized to [0,1])
    %   processed - Processed image (same size as original)
    % Output:
    %   psnr_value - PSNR value in decibels (dB)

    % Ensure the inputs are of type double
    original = double(original);
    processed = double(processed);

    % Calculate Mean Squared Error (MSE)
    mse = mean((original(:) - processed(:)).^2);

    % Determine maximum pixel value for normalized images
    max_pixel_value = 1; % Since images are normalized to [0, 1]

    % Calculate PSNR
    if mse == 0
        psnr_value = Inf; % If MSE is 0, PSNR is infinite
    else
        psnr_value = 10 * log10((max_pixel_value^2) / mse);
    end
end
