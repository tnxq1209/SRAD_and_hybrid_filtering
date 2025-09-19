% THIS IS JUST A TEMPLATE
% Final results to put together in a spread sheet
    analysis = {'PSNR', 'RMSE', 'SSIM', 'Correlation'};
    NOISE = [psnr_noisy, rmse_noisy,ssim_noisy,cp_noisy];
    FILT = [psnr_lee, rmse_lee,ssim_lee,cp_lee];
  

    % Put results into a table
    T = table(analysis', NOISE', FILT'', 'VariableNames', {'ANALYSIS', 'NOISY', 'FILT'});

    % Write to Excel file
    writetable(T, 'Filter_Comparison.xlsx');

    disp('Results saved to Filter_Comparison.xlsx');