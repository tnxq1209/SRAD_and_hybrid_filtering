% Main Script for Speckle Noise Reduction with Adjustable Parameters (2016a)
% Parameters for SRAD
kappa = 10;             % Conduction coefficient
lambda = 0.1;           % Time step size
iterations = 10;        % Number of iterations
epsilon_srad = 1e-6;    % Small constant to prevent division by zero

% Parameters for Wavelet Decomposition
wavelet_type = 'db1';   % Wavelet type for decomposition
wavelet_level = 3;      % Number of decomposition levels

% Parameters for Guided Filter
radius = 5;             % Radius for guided filter (in pixels)
epsilon_guided = 1e-5;  % Regularization parameter for guided filter

% Load the input image (replace with your image file path)
img = imread('C:\Users\Dell\Desktop\Major_Project\h60.jpg'); 
if size(img, 3) == 3
    img = rgb2gray(img);  % Convert to grayscale if RGB
end
img = im2double(img); % Convert to double precision and normalize to [0,1]

% Add speckle noise for demonstration
noise_level = 0.08;
noisy_img = imnoise(img, 'speckle', noise_level);

% Step 1: Apply SRAD Filter
srad_filtered_img = srad_filter(noisy_img, kappa, lambda, iterations, epsilon_srad);

% Step 2: Apply Wavelet Decomposition and Reconstruction
wavelet_filtered_img = wavelet_decomposition(srad_filtered_img, wavelet_type, wavelet_level);

% Step 3: Apply Guided Filter
guided_filtered_img = guided_filter(wavelet_filtered_img, img, radius, epsilon_guided);

% Step 4: Calculate PSNR for each step
psnr_noisy = calculate_psnr(img, noisy_img);
psnr_srad = calculate_psnr(img, srad_filtered_img);
psnr_wavelet = calculate_psnr(img, wavelet_filtered_img);
psnr_guided = calculate_psnr(img, guided_filtered_img);

% Step 5: Calculate RMSE for each step
rmse_noisy = calculate_rmse(img, noisy_img);
rmse_srad = calculate_rmse(img, srad_filtered_img);
rmse_wavelet = calculate_rmse(img, wavelet_filtered_img);
rmse_guided = calculate_rmse(img, guided_filtered_img);

% Display the results
figure;
subplot(1, 2, 1); 
imshow(noisy_img, []); 
title(sprintf('Noisy Image\nPSNR: %.2f dB, RMSE: %.5f', psnr_noisy, rmse_noisy));

subplot(1, 2, 2); 
imshow(guided_filtered_img, []); 
title(sprintf('Hybrid Filtered Image\nPSNR: %.2f dB, RMSE: %.5f', psnr_guided, rmse_guided));


