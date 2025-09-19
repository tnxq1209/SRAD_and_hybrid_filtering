% PSNR values for each filter
psnr_values = struct( ...
    'Lee', [31.5, 28.95, 28.19, 32.92, 30.99, 26.81, 32.81, 31.95, 32.2, 28.35, 28.92, 28.03, 27.44, 27.89, 31.1], ...
    'Kuan', [31.26, 28.39, 27.2, 32.48, 30.08, 25.45, 32.33, 31.16, 31.26, 27.25, 27.63, 26.8, 26.22, 26.69, 30.61], ...
    'Frost', [39.52, 37.38, 36.46, 39.78, 37.62, 35.07, 39.82, 38.04, 38.07, 35.41, 35.48, 35.69, 35.48, 35.25, 38.98], ...
    'SRAD', [38.55, 37, 36.06, 39.01, 36.65, 35, 38.8, 36.87, 37.34, 34.75, 34.62, 35.15, 35.15, 34.79, 38.12], ...
    'Wavelet', [31.48, 28.43, 27.39, 33.18, 30.7, 25.43, 33.1, 31.999, 32.33, 27.52, 28, 27.03, 26.28, 26.85, 30.92], ...
    'Guided', [34.71, 31.39, 30.57, 36.9, 34.38, 28.79, 36.67, 35.44, 35.58, 30.66, 31.46, 30.54, 29.69, 30.19, 34.45], ...
    'Hybrid', [38.18, 37.08, 36.86, 40.3, 37.86, 35.27, 39.96, 38.14, 38.17, 35.86, 35.73, 36.01, 35.86, 35.64, 39.34] ...
);

% Image numbers
images = 1:15;

% Extract filter names
filters = fieldnames(psnr_values);

% Line graph
figure;
hold on;
for i = 1:length(filters)
    plot(images, psnr_values.(filters{i}), '-o', 'LineWidth', 1.5);
end

% Set x-axis ticks to correspond to image numbers (1 to 15)
set(gca, 'xtick', images);

% Customize labels and title
xlabel('Image Number');
ylabel('PSNR (dB)');
title('PSNR Comparison (Line Graph)');
legend(filters, 'Location', 'best');
grid on;
hold off;
