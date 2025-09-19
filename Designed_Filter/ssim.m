% this code has no meaning on its own!!
    % Calculate SSIM for ALL filters
    [ssim_noisy, ~]  = ssim(noisy_image, image);
    [ssim_lee, ~]    = ssim(lee_filtered, image);
    [ssim_kuan, ~]   = ssim(kuan_filtered, image);
    [ssim_frost, ~]  = ssim(frost_filtered, image);

    % Display results
    figure;
    subplot(2,2,1); imshow(noisy_image, []); 
    title(sprintf('Noisy Image\nPSNR: %.2f dB, RMSE: %.5f, SSIM:%.4f', psnr_noisy, rmse_noisy,ssim_noisy));
    
    subplot(2,2,2); imshow(lee_filtered, []); 
    title(sprintf('Lee Filter\nPSNR: %.2f dB, RMSE: %.5f, SSIM:%.4f', psnr_lee, rmse_lee,ssim_lee));

    subplot(2,2,3); imshow(kuan_filtered, []); 
    title(sprintf('Kuan Filter\nPSNR: %.2f dB, RMSE: %.5f, SSIM:%.4f', psnr_kuan, rmse_kuan,ssim_kuan));

    subplot(2,2,4); imshow(frost_filtered, []); 
    title(sprintf('Frost Filter\nPSNR: %.2f dB, RMSE: %.5f, SSIM:%.4f', psnr_frost, rmse_frost,ssim_frost));