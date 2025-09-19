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

    % Calculate PSNR values
    psnr_noisy = calculate_psnr(img, noisy_img);
    psnr_filtered = calculate_psnr(img, filtered_img);

    % Calculate RMSE values
    rmse_noisy = calculate_rmse(img, noisy_img);
    rmse_filtered = calculate_rmse(img, filtered_img);

    % Calculate SSIM values
    [ssim_noisy, ~] = ssim(noisy_img, img);
    [ssim_filtered, ~] = ssim(filtered_img, img);

    %Correlation Parameter
    cp_noisy = corr2(noisy_img, img);
    cp_filt = corr2(filtered_img, img);

    % Display Bar graphs
    filters = {'NOISY', 'SRAD'};
    PSNR = [psnr_noisy,psnr_filtered];
    RMSE = [rmse_noisy,rmse_filtered];
    SSIM = [ssim_noisy,ssim_filtered];
    CORP = [cp_noisy,cp_filt];
    data = [PSNR' RMSE' SSIM' CORP'];   % combine columns
    bar(data);
    set(gca, 'XTickLabel', filters);
    legend({'PSNR','RMSE','SSIM','CORP'});
    ylabel('Value');
    title('Filter Performance Comparison');
    saveas(gcf,'bargraph_SRAD.png');

    % Display the results
    figure;
    subplot(1, 3, 1); imshow(img, []); title('Original Image');
    subplot(1, 3, 2); imshow(noisy_img, []); 
    title(sprintf('Noisy Image\nPSNR: %.2f dB\nRMSE: %.5f\nSSIM: %.4f\nCP: %.4f', psnr_noisy, rmse_noisy,ssim_noisy,cp_noisy));
    subplot(1, 3, 3); imshow(filtered_img, []); 
    title(sprintf('SRAD Filtered Image\nPSNR: %.2f dB\nRMSE: %.5f\nSSIM: %.4f\nCP: %.4f', psnr_filtered, rmse_filtered,ssim_filtered,cp_filt));
    saveas(gcf, 'SRAD_img.png');

    % Final results to put together in a spread sheet
    analysis = {'PSNR', 'RMSE', 'SSIM', 'Correlation'};
    NOISE = [psnr_noisy, rmse_noisy,ssim_noisy,cp_noisy];  
    SRAD = [psnr_filtered, rmse_filtered,ssim_filtered,cp_filt]; 

    % Put results into a table
    T = table(analysis', NOISE', SRAD' ,'VariableNames', {'ANALYSIS', 'NOISE', 'SRAD Filtered'});

    % Write to Excel file
    writetable(T, 'Filter_Comparison_SRAD.xlsx');

    disp('Results saved to Filter_Comparison_SRAD.xlsx');


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

function psnr_value = calculate_psnr(original, processed)
    % Calculate PSNR
    original = double(original);
    processed = double(processed);
    mse = mean((original(:) - processed(:)).^2);
    max_pixel_value = 1; % Since the image is normalized to [0, 1]
    psnr_value = 10 * log10((max_pixel_value^2) / mse);
end

function rmse_value = calculate_rmse(original, processed)
    % Calculate RMSE
    original = double(original);
    processed = double(processed);
    rmse_value = sqrt(mean((original(:) - processed(:)).^2));
end

