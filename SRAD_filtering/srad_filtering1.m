function srad_filtering1()
    % SRAD Filtering on an Image with PSNR and RMSE Calculation
    % Parameters
    kappa = 10;          % Conduction coefficient (adjustable)
    lambda = 0.1;        % Time step size
    iterations = 10;     % Number of iterations
    epsilon = 1e-6;      % Small constant to prevent division by zero (adjustable)

    % Load the image
    img = imread('C:\Users\TNXQ\Desktop\Speckle_noise_reduction_using_various_filters\images\h60.jpg'); % Replace with your image path
    if size(img, 3) == 3
        img = rgb2gray(img);  % Convert to grayscale if the image is RGB
    end
    img = double(img) / 255;  % Normalize to [0, 1]

    % Add speckle noise
    noise_level = 0.08; % Speckle noise variance
    noisy_img = imnoise(img, 'speckle', noise_level);

    % Apply SRAD
    filtered_img = srad_filter(noisy_img, kappa, lambda, iterations, epsilon);

    % Calculate SSIM values
    [ssim_noisy, ~] = ssim(noisy_img, img);
    [ssim_filtered, ~] = ssim(filtered_img, img);

    %Correlation Parameter
    cp_noisy = corr2(noisy_img, img);
    cp_filt = corr2(filtered_img, img);


    % Display the results
    figure;
    subplot(1, 3, 1); imshow(img, []); title('Original Image');
    subplot(1, 3, 2); imshow(noisy_img, []); 
    title(sprintf('Noisy Image\nPSNR: %.2f dB\nRMSE: %.5f\nSSIM: %.4f\nCP: %.4f', ssim_noisy,cp_noisy));
    subplot(1, 3, 3); imshow(filtered_img, []); 
    title(sprintf('SRAD Filtered Image\nPSNR: %.2f dB\nRMSE: %.5f\nSSIM: %.4f\nCP: %.4f',ssim_filtered,cp_filt));
    saveas(gcf, 'SRAD_img.png');

end

function filtered_img = srad_filter(img, kappa, lambda, iterations, epsilon)
    % SRAD Filtering Function with adjustable kappa and epsilon
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
        q0_squared = variance / (mean_squared^2 + epsilon); % Adjusted with epsilon

        c_north = exp(-((north ./ (filtered_img + epsilon)).^2) ./ (kappa * q0_squared + epsilon));
        c_south = exp(-((south ./ (filtered_img + epsilon)).^2) ./ (kappa * q0_squared + epsilon));
        c_west = exp(-((west ./ (filtered_img + epsilon)).^2) ./ (kappa * q0_squared + epsilon));
        c_east = exp(-((east ./ (filtered_img + epsilon)).^2) ./ (kappa * q0_squared + epsilon));

        % Update the image using diffusion equation
        diffusion = lambda * (c_north .* north + c_south .* south + ...
                              c_west .* west + c_east .* east);
        filtered_img = filtered_img + diffusion;
    end
end