function snr_val = calculate_snr(original, processed)
    % calculate_snr - Calculates the Signal-to-Noise Ratio (SNR) 
    % between the original and processed images.
    % Inputs:
    %   processed - Processed (noisy or filtered) image.
    %   original  - Original (reference) image.
    % Output:
    %   snr_val  - Signal-to-Noise Ratio (SNR) in dB.

    % Calculate signal power (mean of original image squared)
    signal_power = mean(original(:).^2);  
    
    % Calculate noise power (mean of squared difference between original and processed)
    noise_power = mean((original(:) - processed(:)).^2);  
    
    % Calculate SNR in dB
    snr_val = 10 * log10(signal_power / noise_power);
end
