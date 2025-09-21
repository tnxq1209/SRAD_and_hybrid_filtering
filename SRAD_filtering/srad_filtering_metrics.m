function srad_filtering1()
    % SRAD Filtering on an Image with PSNR and RMSE Calculation
    % Parameters
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
        img = imread(fpath);% reads the input images.
    
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
        ssim_value = ssim(filtered_img, img);
    
        %Correlation Parameter
        cp_value = corr2(filtered_img, img);

        % === Store in Tables ===
        newRowSSIM = {fname, ssim_value};
        SSIM_Table = [SSIM_Table; newRowSSIM]; %#ok<AGROW>

        newRowCorr = {fname, cp_value};
        Corr_Table = [Corr_Table; newRowCorr]; %#ok<AGROW>

    end
    
    % === 5. Add Column Names ===
    SSIM_Table.Properties.VariableNames = { ...
        'ImageName','SSIM_SRAD'};
    Corr_Table.Properties.VariableNames = { ...
        'ImageName','Corr_SRAD'};

    % === 6. Save to Excel ===
    writetable(SSIM_Table, fullfile(outFolder, 'SSIM_Results.xlsx'));
    writetable(Corr_Table, fullfile(outFolder, 'Correlation_Results.xlsx'));

    disp('All images processed and results saved!');


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