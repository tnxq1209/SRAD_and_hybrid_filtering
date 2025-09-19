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
    saveas(gcf,'bargraph.png');
