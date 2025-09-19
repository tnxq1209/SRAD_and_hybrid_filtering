function wavelet_filtering1()
    % Step 1: Read and preprocess the image
    img = imread('C:\Users\TNXQ\Desktop\Speckle_noise_reduction_using_various_filters-main\Speckle_noise_reduction_using_various_filters\Speckle_noise_reduction_using_various_filters\images\h48.jpg');  % Replace with your image path
    if size(img, 3) == 3
        img = rgb2gray(img);  % Convert to grayscale if the image is RGB
    end
    img = double(img) / 255;  % Normalize the image to [0, 1]

    % Step 2: Add Speckle Noise (0.08 variance noise level)
    noise_img = imnoise(img, 'speckle', 0.08);  % Adding speckle noise with variance

    % Step 3: Perform 2D Wavelet Decomposition
    level = 1;  % Decomposition level (use level 1 for simplicity)
    waveletType = 'db4';  % Daubechies 4 wavelet
    [c, s] = wavedec2(noise_img, level, waveletType);  % Decompose the noisy image

    % Step 4: Extract Approximation and Detail Coefficients
    LL = appcoef2(c, s, waveletType, level);  % Approximation coefficients (LL)
    [LH, HL, HH] = detcoef2('all', c, s, level);  % Detail coefficients (LH, HL, HH)

    % Step 5: Apply Thresholding to the Detail Coefficients (Optional)
    threshold = 0.1;  % Define a threshold value for filtering
    LH = LH .* (abs(LH) > threshold);  % Thresholding horizontal details
    HL = HL .* (abs(HL) > threshold);  % Thresholding vertical details
    HH = HH .* (abs(HH) > threshold);  % Thresholding diagonal details

    % Step 6: Reconstruct the Wavelet Coefficients
    c_filtered = c;  % Initialize with the original coefficients
    c_filtered(1:numel(LL)) = LL(:);  % Update the approximation coefficients (LL)

    % Update the detail coefficients in the correct order
    idx = numel(LL) + 1;
    c_filtered(idx:idx + numel(LH) - 1) = LH(:);
    idx = idx + numel(LH);
    c_filtered(idx:idx + numel(HL) - 1) = HL(:);
    idx = idx + numel(HL);
    c_filtered(idx:idx + numel(HH) - 1) = HH(:);

    % Step 7: Reconstruct the Image from the Filtered Coefficients
    filtered_img = waverec2(c_filtered, s, waveletType);  % Reconstruct the filtered image

    % Step 8: Calculate PSNR and RMSE Values using the custom functions
    psnr_noisy = calculate_psnr(noise_img, img);  % PSNR for noisy image
    rmse_noisy = calculate_rmse(noise_img, img);  % RMSE for noisy image
    psnr_filtered = calculate_psnr(filtered_img, img);  % PSNR for filtered image
    rmse_filtered = calculate_rmse(filtered_img, img);  % RMSE for filtered image

    % Calculating SSIM
    [ssim_noisy, ~]  = ssim(noise_img, img);
    [ssim_filterd, ~]  = ssim(filtered_img, img);

    %Correlation parameter
    cp_noisy = corr2(noise_img, img);
    cp_filt = corr2(filtered_img, img);

    % Step 9: Display the Results
    figure;

    % Original Image
    subplot(1, 3, 1);
    imshow(img, []);
    title('Original Image');

    % Noisy Image
    subplot(1, 3, 2);
    imshow(noise_img, []);
    title(sprintf('Noisy Image\nPSNR: %.2f dB\nRMSE: %.5f\nSSIM: %.4f\nCP: %.4f', psnr_noisy, rmse_noisy,ssim_noisy,cp_noisy));

    % Reconstructed (Filtered) Image
    subplot(1, 3, 3);
    imshow(filtered_img, []);
    title(sprintf('Filtered Image\nPSNR: %.2f dB\nRMSE: %.5f\nSSIM: %.4f\nCP: %.4f', psnr_filtered, rmse_filtered,ssim_filterd,cp_filt));
end

% Custom PSNR Function
function psnr_val = calculate_psnr(processed, original)
    % calculate_psnr - Calculates the Peak Signal-to-Noise Ratio (PSNR) 
    % between the original and processed images.
    %
    % Inputs:
    %   processed - Processed (noisy or filtered) image.
    %   original  - Original (reference) image.
    %
    % Output:
    %   psnr_val  - Peak Signal-to-Noise Ratio (PSNR) in dB.

    mse = mean((processed(:) - original(:)).^2);  % Mean Squared Error
    max_val = 1;  % Maximum possible pixel value (for normalized images, it's 1)
    
    % Calculate PSNR
    psnr_val = 10 * log10(max_val^2 / mse);
end

% Custom RMSE Function
function rmse_val = calculate_rmse(processed, original)
    % calculate_rmse - Calculates the Root Mean Squared Error (RMSE)
    % between the original and processed images.
    %
    % Inputs:
    %   processed - Processed (noisy or filtered) image.
    %   original  - Original (reference) image.
    %
    % Output:
    %   rmse_val  - Root Mean Squared Error (RMSE).

    mse = mean((processed(:) - original(:)).^2);  % Mean Squared Error
    rmse_val = sqrt(mse);  % RMSE is the square root of MSE
end
