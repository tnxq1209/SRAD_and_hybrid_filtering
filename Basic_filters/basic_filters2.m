function basic_filters1()
    % Speckle Noise Reduction using Lee, Kuan, and Frost Filters

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
        
        % Calculate SSIM for ALL filters
        ssim_noisy  = ssim( noisy_image, image);
        ssim_lee    = ssim( lee_filtered, image);
        ssim_kuan   = ssim( kuan_filtered, image);
        ssim_frost = ssim(frost_filtered, image);
    
        % correlation parameter
        cp_noisy = corr2(noisy_image, image);
        cp_lee = corr2(lee_filtered, image);
        cp_kuan = corr2(kuan_filtered, image);
        cp_frost = corr2(frost_filtered, image);
   
        % === Store in Tables ===
        newRowSSIM = {fname, ssim_noisy, ssim_lee, ssim_kuan, ssim_frost};
        SSIM_Table = [SSIM_Table; newRowSSIM]; %#ok<AGROW>

        newRowCorr = {fname, cp_noisy, cp_lee, cp_kuan, cp_frost};
        Corr_Table = [Corr_Table; newRowCorr]; %#ok<AGROW>

    end
    
    % === 5. Add Column Names ===
    SSIM_Table.Properties.VariableNames = { ...
        'ImageName','SSIM_Noisy','SSIM_Lee','SSIM_Kuan','SSIM_Frost'};
    Corr_Table.Properties.VariableNames = { ...
        'ImageName','Corr_Noisy','Corr_Lee','Corr_Kuan','Corr_Frost'};

    % === 6. Save to Excel ===
    writetable(SSIM_Table, fullfile(outFolder, 'SSIM_Results.xlsx'));
    writetable(Corr_Table, fullfile(outFolder, 'Correlation_Results.xlsx'));

    disp('All images processed and results saved!');

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