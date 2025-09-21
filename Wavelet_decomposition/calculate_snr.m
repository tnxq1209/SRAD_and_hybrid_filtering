% Function to Calculate SNR
function snr_value = calculate_snr(original, processed)
    % Ensure images are double for calculations
    original = double(original);
    processed = double(processed);

    % Calculate signal power (mean squared value of the original image)
    signal_power = mean(original(:).^2);

    % Calculate noise power
    noise_power = mean((original(:) - processed(:)).^2);

    % Compute SNR in dB
    snr_value = 10 * log10(signal_power / noise_power);
end

