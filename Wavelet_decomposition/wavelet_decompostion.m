% Step 1: Read and preprocess the image
img = imread('C:\Users\TNXQ\Desktop\Speckle_noise_reduction_using_various_filters-main\Speckle_noise_reduction_using_various_filters\Speckle_noise_reduction_using_various_filters\images\h96.jpg');  % Replace with your image path
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

% Step 6: Reconstruct the Image from the Filtered Coefficients
c_filtered = [LL(:); LH(:); HL(:); HH(:)];  % Recombine the coefficients
filtered_img = waverec2(c_filtered, s, waveletType);  % Reconstruct the filtered image

% Step 7: Calculate SNR Values
snr_noisy = calculate_snr(img, noise_img);  % SNR for noisy image
snr_filtered = calculate_snr(img, filtered_img);  % SNR for filtered image

% Calculating SSIM
[ssim_noisy, ~]  = ssim(noise_img, img);
[ssim_filterd, ~]  = ssim(filtered_img, img);

%Correlation parameter
cp_noisy = corr2(noise_img, img);
cp_filt = corr2(filtered_img, img);

% Step 8: Display the Results
figure;

% Original Image
subplot(1, 3, 1);
imshow(img, []);
title('Original Image');

% Noisy Image
subplot(1, 3, 2);
imshow(noise_img, []);
title(sprintf('Noisy Image\nSNR: %.2f dB\nSSIM: %.4f\nCP: %.4f', snr_noisy,ssim_noisy,cp_noisy));

% Reconstructed (Filtered) Image
subplot(1, 3, 3);
imshow(filtered_img, []);
title(sprintf('Filtered Image\nSNR: %.2f dB\nSSIM: %.4f\nCP:%.4f', snr_filtered,ssim_filterd,cp_filt));

