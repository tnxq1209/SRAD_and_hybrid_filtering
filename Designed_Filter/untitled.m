% === 1. Define Input & Output Paths ===
inFolder  = 'C:\Users\TNXQ\Desktop\Speckle_noise_reduction_using_various_filters\img_2';   % input images folder
outFolder = 'C:\Users\TNXQ\Desktop\Speckle_noise_reduction_using_various_filters\Metrices_table';       % output results folder

% Create output folder if not exists
if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

% === 2. Get List of Images ===
imgFiles = dir(fullfile(inFolder, '*.jpg'));   % change to *.jpg if needed

% === 3. Prepare Results Table ===
Results = table;

% === 4. Loop Through Images ===
for k = 1:length(imgFiles)
    % Read image
    fname = imgFiles(k).name;
    fpath = fullfile(inFolder, fname);
    origImg = imread(fpath);

    % Convert to grayscale if needed
    if size(origImg,3) == 3
        origImg = rgb2gray(origImg);
    end

    % Simulate noisy image (if you already have noisy images, skip this step)
    noisyImg = imnoise(origImg, 'speckle', 0.04);

    % === Save Filtered Images (Optional) ===
    imwrite(filtered_lee,   fullfile(outFolder, [fname '_Lee.png']));
    imwrite(filtered_kuan,  fullfile(outFolder, [fname '_Kuan.png']));
    imwrite(filtered_frost, fullfile(outFolder, [fname '_Frost.png']));

    % === Store Results in Table ===
    newRow = {
        fname, ...
        ssim_noisy, corr_noisy, ...
        ssim_lee, corr_lee, ...
        ssim_kuan, corr_kuan, ...
        ssim_frost, corr_frost
    };

    Results = [Results; newRow]; %#ok<AGROW>
end

% === 5. Add Column Names ===
Results.Properties.VariableNames = { ...
    'ImageName', ...
    'SSIM_Noisy','Corr_Noisy', ...
    'SSIM_Lee','Corr_Lee', ...
    'SSIM_Kuan','Corr_Kuan', ...
    'SSIM_Frost','Corr_Frost' ...
};

% === 6. Save to Excel ===
writetable(Results, fullfile(outFolder, 'Filter_Results.xlsx'));

disp('All images processed and results saved!');


