function guided_filter_with_metrics()
    % Guided Filter with PSNR and RMSE Calculation for Speckle Noise

    clc; clear; close all;

    % Load an image
    image = imread('C:\Users\TNXQ\Desktop\Speckle_noise_reduction_using_various_filters\images\h60.jpg'); % Replace with your image file
    if size(image, 3) == 3
        image = rgb2gray(image); % Convert to grayscale if it's a color image
    end

    % Normalize the image to range [0, 1]
    image = double(image) / 255;

    % Add speckle noise with variance 0.08
    noise_variance = 0.08;
    noisy_image = imnoise(image, 'speckle', noise_variance);

    % Parameters for the guided filter
    radius = 5;          % Radius of the window
    epsilon = 0.02;    % Regularization parameter

    % Apply guided filter
    filtered_image = guided_filter(noisy_image, noisy_image, radius, epsilon);

    % Calculate PSNR and RMSE and SSIM and CP
    psnr_noisy = calculate_psnr(noisy_image, filtered_image);
    rmse_noisy = calculate_rmse(noisy_image, filtered_image);
    [ssim_noisy, ~]  = ssim(noisy_image, filtered_image);
    cp_noisy = corr2(noisy_image, filtered_image);

    psnr_value = calculate_psnr(image, filtered_image);
    rmse_value = calculate_rmse(image, filtered_image);
    [ssim_value, ~]  = ssim(image, filtered_image);
    cp_value = corr2(image, filtered_image);


    % Display Bar graphs
    filters = {'NOISY', 'FILTERED'};
    PSNR = [psnr_noisy,psnr_value];
    RMSE = [psnr_noisy,rmse_value];
    SSIM = [ssim_noisy,ssim_value];
    CORP = [cp_noisy,cp_value];
    data = [PSNR' RMSE' SSIM' CORP'];   % combine columns
    bar(data);
    set(gca, 'XTickLabel', filters);
    legend({'PSNR','RMSE','SSIM','CORP'});
    ylabel('Value');
    title('Filter Performance Comparison');
    saveas(gcf,'bargraph_guided.png');

    % Display results
    figure;
    subplot(1, 2, 1); imshow(noisy_image, []); title(sprintf('Noisy Image\n(Speckle Variance: %.2f)', noise_variance));
    subplot(1, 2, 2); imshow(filtered_image, []);
    title(sprintf('Guided Filter\nPSNR: %.2f dB, RMSE: %.5f, SSIM: %.4f, CP:%.4f', psnr_value, rmse_value,ssim_noisy,cp_noisy));
    % Final results to put together in a spread sheet
    analysis = {'PSNR', 'RMSE', 'SSIM', 'Correlation'};
    NOISE = [psnr_value, rmse_value,ssim_noisy,cp_noisy];  
    saveas(gcf,'guided_filtered_img.png');


    % Put results into a table
    T = table(analysis', NOISE', 'VariableNames', {'ANALYSIS', 'GUIDED FILTER'});

    % Write to Excel file
    writetable(T, 'Filter_Comparison_guided.xlsx');

    disp('Results saved to Filter_Comparison_guided.xlsx');

end

function guided_img = guided_filter(I, P, radius, epsilon)
    % I: Guidance image (input image)
    % P: Filtering image (input image)
    % radius: Window radius
    % epsilon: Regularization parameter

    [rows, cols] = size(I);
    N = box_filter(ones(rows, cols), radius);

    mean_I = box_filter(I, radius) ./ N;
    mean_P = box_filter(P, radius) ./ N;
    corr_I = box_filter(I .* I, radius) ./ N;
    corr_IP = box_filter(I .* P, radius) ./ N;

    var_I = corr_I - mean_I .* mean_I;
    cov_IP = corr_IP - mean_I .* mean_P;

    a = cov_IP ./ (var_I + epsilon);
    b = mean_P - a .* mean_I;

    mean_a = box_filter(a, radius) ./ N;
    mean_b = box_filter(b, radius) ./ N;

    guided_img = mean_a .* I + mean_b;
end

function result = box_filter(img, radius)
    % Fast box filter for mean calculation
    kernel = ones(2*radius + 1, 2*radius + 1);
    result = conv2(img, kernel, 'same');
end

function psnr_value = calculate_psnr(original, processed)
    % PSNR Calculation
    mse = mean((original(:) - processed(:)).^2);
    max_pixel_value = 1; % Since images are normalized to [0, 1]
    psnr_value = 10 * log10((max_pixel_value^2) / mse);
end

function rmse_value = calculate_rmse(original, processed)
    % RMSE Calculation
    mse = mean((original(:) - processed(:)).^2);
    rmse_value = sqrt(mse);
end
