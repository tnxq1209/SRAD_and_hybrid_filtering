function wavelet_filtering1()
    % Step 1: Read and preprocess the image
    img = imread('C:\Users\TNXQ\Desktop\Speckle_noise_reduction_using_various_filters\images\h48.jpg');  % Replace with your image path
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
    title(sprintf('Filtered Image\nSSIM: %.4f\nCP: %.4f', ssim_filterd,cp_filt));
    saveas(gcf,'wavelet_filtering.png');
end