function basic_filters()
    % Speckle Noise Reduction using Lee, Kuan, and Frost Filters

    clc; clear; close all;

    % Load an image with speckle noise
    image = imread('C:\Users\TNXQ\Desktop\Speckle_noise_reduction_using_various_filters-main\Speckle_noise_reduction_using_various_filters\Speckle_noise_reduction_using_various_filterss\h96.jpg'); % Replace with your image file
    if size(image, 3) == 3
        image = rgb2gray(image); % Convert to grayscale if it's a color image
    end

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

    % Calculate SNR values
    snr_noisy = calculate_snr(image, noisy_image);
    snr_lee = calculate_snr(image, lee_filtered);
    snr_kuan = calculate_snr(image, kuan_filtered);
    snr_frost = calculate_snr(image, frost_filtered);

    % Display results
    figure;
    subplot(2,2,1); imshow(noisy_image, []); title(sprintf('Noisy Image\nSNR: %.2f dB', snr_noisy));
    subplot(2,2,2); imshow(lee_filtered, []); title(sprintf('Lee Filter\nSNR: %.2f dB', snr_lee));
    subplot(2,2,3); imshow(kuan_filtered, []); title(sprintf('Kuan Filter\nSNR: %.2f dB', snr_kuan));
    subplot(2,2,4); imshow(frost_filtered, []); title(sprintf('Frost Filter\nSNR: %.2f dB', snr_frost));
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

% ----- SNR Calculation Function -----
function snr_value = calculate_snr(original, processed)
    original = double(original);
    processed = double(processed);
    signal_power = sum(original(:).^2);
    noise_power = sum((original(:) - processed(:)).^2);
    snr_value = 10 * log10(signal_power / noise_power);
end
