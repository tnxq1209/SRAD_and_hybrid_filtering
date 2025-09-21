function guided_filter_with_metrics()
    % Guided Filter with PSNR and RMSE Calculation for Speckle Noise

    clc; clear; close all;

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
    
        % Calculate SSIM and CP
        [ssim_value, ~]  = ssim(image, filtered_image);
        cp_value = corr2(image, filtered_image);
       
        % === Store in Tables ===
        newRowSSIM = {fname, ssim_value};
        SSIM_Table = [SSIM_Table; newRowSSIM]; %#ok<AGROW>

        newRowCorr = {fname, cp_value};
        Corr_Table = [Corr_Table; newRowCorr]; %#ok<AGROW>

    end
    
    % === 5. Add Column Names ===
    SSIM_Table.Properties.VariableNames = { ...
        'ImageName','SSIM_Guided'};
    Corr_Table.Properties.VariableNames = { ...
        'ImageName','Corr_Guided'};

    % === 6. Save to Excel ===
    writetable(SSIM_Table, fullfile(outFolder, 'SSIM_Results.xlsx'));
    writetable(Corr_Table, fullfile(outFolder, 'Correlation_Results.xlsx'));

    disp('All images processed and results saved!');

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