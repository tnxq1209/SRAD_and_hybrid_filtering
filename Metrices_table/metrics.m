function final_filtering_all()
    % Speckle Noise Reduction using Lee, Kuan, and Frost Filters

    clc; clear; close all;

    % Parameters for srad.
    kappa = 10;          % Conduction coefficient (adjustable)
    lambda = 0.1;        % Time step size
    iterations = 10;     % Number of iterations
    epsilon = 1e-6;      % Small constant to prevent division by zero (adjustable)

    % === 1. Define Input & Output Paths ===
    inFolder  = 'C:\Users\TNXQ\Desktop\Speckle_noise_reduction_using_various_filters\images\image_inuse';   % input images folder(choose your own image folder)
    outFolder = 'C:\Users\TNXQ\Desktop\Speckle_noise_reduction_using_various_filters\Metrices_table';       % output results folder(choose your own output folder)

    % Create output folder if not exists
    if ~exist(outFolder, 'dir')
        mkdir(outFolder);
    end

    % === 2. Get List of Images ===
    imgFiles = dir(fullfile(inFolder, '*.jpg'));   % change to *.png if needed

    % === 3. Prepare Results Table ===
    SSIM_Table = table; % XLSX table where result of all images with SSIM values.
    Corr_Table = table; % XLSX table where result of all images with Corelation parameter values.


    % === 4. Loop Through Images ===
    for k = 1:length(imgFiles)
        % Read image
        fname = imgFiles(k).name;% extract file names from the list of images.
        fpath = fullfile(inFolder, fname);% combine input folder and image name to get full path where image is stored.
        image = imread(fpath);% reads the input images.
    
        % Convert to grayscale if needed
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

        % Parameters for the guided filter
        radius = 5;          % Radius of the window
        epsilon = 0.02;    % Regularization parameter
    
        % Apply guided filter
        guided_image = guided_filter(noisy_image, noisy_image, radius, epsilon);

        % WAVELET: Perform 2D Wavelet Decomposition
        level = 1;  % Decomposition level (use level 1 for simplicity)
        waveletType = 'db4';  % Daubechies 4 wavelet
        [c, s] = wavedec2(noisy_image, level, waveletType);  % Decompose the noisy image
    
        % WAVELET: Extract Approximation and Detail Coefficients
        LL = appcoef2(c, s, waveletType, level);  % Approximation coefficients (LL)
        [LH, HL, HH] = detcoef2('all', c, s, level);  % Detail coefficients (LH, HL, HH)
    
        % WAVELET: Apply Thresholding to the Detail Coefficients (Optional)
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
    
        % WAVELET: Reconstruct the Image from the Filtered Coefficients
        wavelet_img = waverec2(c_filtered, s, waveletType);  % Reconstruct the filtered image

        % Apply SRAD
        srad_img = srad_filter(noisy_image, kappa, lambda, iterations, epsilon);
        
        % Calculate SSIM for ALL filters
        ssim_noisy  = ssim( noisy_image, image);
        ssim_lee    = ssim( lee_filtered, image);
        ssim_kuan   = ssim( kuan_filtered, image);
        ssim_frost = ssim(frost_filtered, image);
        ssim_gd = ssim(guided_image ,image);
        ssim_wav  = ssim(wavelet_img, image);
        ssim_srad = ssim(srad_img, image);

        % correlation parameter
        cp_noisy = corr2(noisy_image, image);
        cp_lee = corr2(lee_filtered, image);
        cp_kuan = corr2(kuan_filtered, image);
        cp_frost = corr2(frost_filtered, image);
        cp_gd = corr2( guided_image,image);
        cp_wav = corr2(wavelet_img, image);
        cp_srad = corr2(srad_img, image);

        % === Store in Tables ===
        newRowSSIM = {fname, ssim_noisy, ssim_lee, ssim_kuan, ssim_frost,ssim_wav,ssim_gd,ssim_srad};
        SSIM_Table = [SSIM_Table; newRowSSIM]; %#ok<AGROW>

        newRowCorr = {fname, cp_noisy, cp_lee, cp_kuan, cp_frost,cp_wav,cp_gd,cp_srad};
        Corr_Table = [Corr_Table; newRowCorr]; %#ok<AGROW>

    end
    
    % === 5. Add Column Names ===
    SSIM_Table.Properties.VariableNames = { ...
        'ImageName','SSIM_NOISY','SSIM_LEE','SSIM_KUAN','SSIM_FROST','SSIM_Wavelet','SSIM_Guided','SSIM_SRAD'};
    Corr_Table.Properties.VariableNames = { ...
        'ImageName','Corr_NOISY','Corr_LEE','Corr_KUAN','Corr_FROST','Corr_Wavelet','Corr_Guided','Corr_SRAD'};

    % === 6. Save to Excel ===
    writetable(SSIM_Table, fullfile(outFolder, 'SSIM_Results.xlsx'));
    writetable(Corr_Table, fullfile(outFolder, 'Correlation_Results.xlsx'));

    disp('All images processed and results saved!');


end


%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxx
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