% Load the input image (replace with your image file)
img = imread('C:\Users\TNXQ\Desktop\Speckle_noise_reduction_using_various_filters-main\Speckle_noise_reduction_using_various_filters\Speckle_noise_reduction_using_various_filters\h96.jpg'); 
img = im2double(img); % Convert to double precision

% Add speckle noise to the image
noisy_img = imnoise(img, 'speckle', 0.08);

% Apply Guided Filter
radius = 5;      % Radius of the guided filter
epsilon = 0.01;  % Regularization parameter
guided_img = guided_filter(noisy_img, radius, epsilon);

% Calculate SNR values
snr_original = snr(img, img - img);           % SNR of the original image (theoretical max)
snr_noisy = snr(img, noisy_img - img);        % SNR of the noisy image
snr_guided = snr(img, guided_img - img);      % SNR of the guided-filtered image

% Display the original, noisy, and guided-filtered images with SNR values
figure;

subplot(1, 3, 1);
imshow(img, []);
title(['Original Image']);

subplot(1, 3, 2);
imshow(noisy_img, []);
title(['Noisy Image (SNR: ', num2str(snr_noisy, '%.2f'), ' dB)']);

subplot(1, 3, 3);
imshow(guided_img, []);
title(['Guided Filtered Image (SNR: ', num2str(snr_guided, '%.2f'), ' dB)']);
