function wavelet_filtering1()
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
        img = double(img) / 255;  % Normalize the image to [0, 1]
    
        % Step 2: Add Speckle Noise (0.08 variance noise level)
        noise_img = imnoise(img, 'speckle', 0.08);  % Adding speckle noise with variance
    
        % Step 3: Perform 2D Wavelet Decomposition
        level = 1;  % Decomposition level (use level 1 for simplicity)
        waveletType = 'db4';  % Daubechies 4 wavelet
        [c, s] = wavedec2(noise_img, level, waveletType);  % Decompose the noisy image
    
        % Step 4: Extract Approximation and Detail Coefficients
        LL = appcoef2(c, s, waveletType, level);  % Approximation coefficients (LL)
        [LH, HL, HH] = detcoef2('all', c, s, level);  % Detail coefficients (LH, HL, HH)
    
        % Step 5: Apply Thresholding to the Detail Coefficients (Optional)
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
    
        % Step 7: Reconstruct the Image from the Filtered Coefficients
        filtered_img = waverec2(c_filtered, s, waveletType);  % Reconstruct the filtered image
    
       
        % Calculating SSIM
        
        ssim_value  = ssim(filtered_img, img);
    
        %Correlation parameter
       
        cp_value = corr2(filtered_img, img);
    
        % === Store in Tables ===
        newRowSSIM = {fname, ssim_value};
        SSIM_Table = [SSIM_Table; newRowSSIM]; %#ok<AGROW>
    
        newRowCorr = {fname, cp_value};
        Corr_Table = [Corr_Table; newRowCorr]; %#ok<AGROW>

    end
    
    % === 5. Add Column Names ===
    SSIM_Table.Properties.VariableNames = { ...
        'ImageName','SSIM_Wavelet'};
    Corr_Table.Properties.VariableNames = { ...
        'ImageName','Corr_Wavelet'};

    % === 6. Save to Excel ===
    writetable(SSIM_Table, fullfile(outFolder, 'SSIM_Results.xlsx'));
    writetable(Corr_Table, fullfile(outFolder, 'Correlation_Results.xlsx'));

    disp('All images processed and results saved!');


end