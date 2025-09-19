function basic_filters1()
    % Speckle Noise Reduction using Lee, Kuan, and Frost Filters

    clc; clear; close all;

    % Load an image with speckle noise
    image = imread('C:\Users\TNXQ\Desktop\Speckle_noise_reduction_using_various_filters-main\Speckle_noise_reduction_using_various_filters\Speckle_noise_reduction_using_various_filters\h60.jpg'); % Replace with your image file
    
    if size(image, 3) == 3
        image = rgb2gray(image); % Convert to grayscale if it's a color image
    end

    % Normalize image to [0, 1]
    image = double(image) / 255;

    % Add speckle noise with 0.08 variance (equivalent to 0.08 dB noise)
    noise_variance = 0.08;
    noisy_image = imnoise(image, 'speckle', noise_variance);

    % Parameters for the filters
    window_size = 5;      % Size of the local filtering window
    damping_factor = 1;   % Damping factor for the Frost filter

    % Step 1: Apply the Lee Filter
    lee_filtered = leeFilter(noisy_image, window_size);

    % Step 2: Apply the Kuan Filter
    kuan_filtered = kuanFilter(noisy_image, window_size);

    % Step 3: Apply the Frost Filter
    frost_filtered = frostFilter(noisy_image, window_size, damping_factor);

    % Calculate PSNR and RMSE values
    [psnr_noisy, rmse_noisy] = calculate_metrics(image, noisy_image);
    [psnr_lee, rmse_lee] = calculate_metrics(image, lee_filtered);
    [psnr_kuan, rmse_kuan] = calculate_metrics(image, kuan_filtered);
    [psnr_frost, rmse_frost] = calculate_metrics(image, frost_filtered);

    % Display results
    figure;
    subplot(2,2,1); imshow(noisy_image, []); 
    title(sprintf('Noisy Image\nPSNR: %.2f dB, RMSE: %.5f', psnr_noisy, rmse_noisy));
    
    subplot(2,2,2); imshow(lee_filtered, []); 
    title(sprintf('Lee Filter\nPSNR: %.2f dB, RMSE: %.5f', psnr_lee, rmse_lee));

    subplot(2,2,3); imshow(kuan_filtered, []); 
    title(sprintf('Kuan Filter\nPSNR: %.2f dB, RMSE: %.5f', psnr_kuan, rmse_kuan));

    subplot(2,2,4); imshow(frost_filtered, []); 
    title(sprintf('Frost Filter\nPSNR: %.2f dB, RMSE: %.5f', psnr_frost, rmse_frost));
end

% ----- Lee Filter Function -----
function filtered_image = leeFilter(img, window_size)
    img = double(img);
    meanLocal = filter2(ones(window_size) / window_size^2, img, 'same');
    varLocal = filter2(ones(window_size) / window_size^2, img.^2, 'same') - meanLocal.^2;
    
    % Estimate noise variance (median of local variances)
    noiseVariance = mean(varLocal(:));
    
    % Lee filter formula
    K = varLocal ./ (varLocal + noiseVariance);
    filtered_image = meanLocal + K .* (img - meanLocal);
end

% ----- Kuan Filter Function -----
function filtered_image = kuanFilter(img, window_size)
    img = double(img);
    meanLocal = filter2(ones(window_size) / window_size^2, img, 'same');
    varLocal = filter2(ones(window_size) / window_size^2, img.^2, 'same') - meanLocal.^2;
    ratio = varLocal ./ (meanLocal .^ 2 + eps);
    K = (1 ./ (1 + ratio));
    filtered_image = meanLocal .* (1 - K) + img .* K;
end

% ----- Frost Filter Function -----
function filtered_image = frostFilter(img, window_size, damping_factor)
    img = double(img);
    [rows, cols] = size(img);
    paddedImg = padarray(img, [floor(window_size/2), floor(window_size/2)], 'symmetric');
    filtered_image = zeros(rows, cols);
    
    for i = 1:rows
        for j = 1:cols
            localWindow = paddedImg(i:i+window_size-1, j:j+window_size-1);
            localMean = mean(localWindow(:));
            coef = exp(-damping_factor * abs(localWindow - localMean) ./ (localMean + eps));
            weights = coef ./ sum(coef(:) + eps);
            filtered_image(i, j) = sum(weights(:) .* localWindow(:));
        end
    end
end

% ----- PSNR and RMSE Calculation Function -----
function [psnr_value, rmse_value] = calculate_metrics(original, processed)
    original = double(original);
    processed = double(processed);

    % Normalize images to [0, 1] if they are not already
    if max(original(:)) > 1
        original = original / 255;
    end
    if max(processed(:)) > 1
        processed = processed / 255;
    end

    % Calculate MSE and RMSE
    mse = mean((original(:) - processed(:)).^2);
    rmse_value = sqrt(mse);

    % Calculate PSNR
    max_pixel_value = 1;  % Since images are normalized to [0, 1]
    psnr_value = 10 * log10((max_pixel_value^2) / mse);
end
