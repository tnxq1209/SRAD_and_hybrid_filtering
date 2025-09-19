function srad_filtering()
    % SRAD Filtering on an Image with SNR Calculation
    % Parameters
    kappa = 20;           % Conduction coefficient
    lambda = 0.1;         % Time step size
    iterations = 20;      % Number of iterations

    % Load the image
    img = imread('C:\Users\Dell\Desktop\Major_Project\h96.jpg'); % Replace with your image path
    if size(img, 3) == 3
        img = rgb2gray(img);  % Convert to grayscale if the image is RGB
    end
    img = double(img) / 255;  % Normalize to [0, 1]

    % Add speckle noise
    noise_level = 0.08; % Speckle noise variance
    noisy_img = imnoise(img, 'speckle', noise_level);

    % Apply SRAD
    filtered_img = srad_filter(noisy_img, kappa, lambda, iterations);

    % Calculate SNR values
    snr_noisy = calculate_snr(img, noisy_img);
    snr_filtered = calculate_snr(img, filtered_img);

    % Display the results
    figure;
    subplot(1, 3, 1); imshow(img, []); title('Original Image');
    subplot(1, 3, 2); imshow(noisy_img, []); 
    title(sprintf('Noisy Image\nSNR: %.2f dB', snr_noisy));
    subplot(1, 3, 3); imshow(filtered_img, []); 
    title(sprintf('SRAD Filtered Image\nSNR: %.2f dB', snr_filtered));
end

function filtered_img = srad_filter(img, kappa, lambda, iterations)
    % SRAD Filtering Function
    [rows, cols] = size(img);
    filtered_img = img; % Initialize filtered image

    for t = 1:iterations
        % Compute gradients
        north = [filtered_img(1, :); filtered_img(1:end-1, :)] - filtered_img;
        south = [filtered_img(2:end, :); filtered_img(end, :)] - filtered_img;
        west = [filtered_img(:, 1), filtered_img(:, 1:end-1)] - filtered_img;
        east = [filtered_img(:, 2:end), filtered_img(:, end)] - filtered_img;

        % Calculate diffusion coefficients
        mean_squared = mean(filtered_img(:).^2);
        variance = mean((filtered_img(:) - mean_squared).^2);
        q0_squared = variance / mean_squared^2;

        c_north = exp(-((north ./ (filtered_img + eps)).^2) ./ (q0_squared + eps));
        c_south = exp(-((south ./ (filtered_img + eps)).^2) ./ (q0_squared + eps));
        c_west = exp(-((west ./ (filtered_img + eps)).^2) ./ (q0_squared + eps));
        c_east = exp(-((east ./ (filtered_img + eps)).^2) ./ (q0_squared + eps));

        % Update the image
        diffusion = lambda * (c_north .* north + c_south .* south + ...
                              c_west .* west + c_east .* east);
        filtered_img = filtered_img + diffusion;
    end
end

function snr_value = calculate_snr(original, processed)
    % Calculate SNR
    original = double(original);
    processed = double(processed);
    signal_power = sum(original(:).^2);
    noise_power = sum((original(:) - processed(:)).^2);
    snr_value = 10 * log10(signal_power / noise_power);
end
