function filtered_img = srad_filter(img, kappa, lambda, iterations, epsilon)
    [rows, cols] = size(img);
    filtered_img = img; % Initialize filtered image

    for t = 1:iterations
        north = [filtered_img(1, :); filtered_img(1:end-1, :)] - filtered_img;
        south = [filtered_img(2:end, :); filtered_img(end, :)] - filtered_img;
        west = [filtered_img(:, 1), filtered_img(:, 1:end-1)] - filtered_img;
        east = [filtered_img(:, 2:end), filtered_img(:, end)] - filtered_img;

        mean_squared = mean(filtered_img(:).^2);
        variance = mean((filtered_img(:) - mean_squared).^2);
        q0_squared = variance / (mean_squared^2 + epsilon);

        c_north = exp(-((north ./ (filtered_img + epsilon)).^2) ./ (kappa * q0_squared + epsilon));
        c_south = exp(-((south ./ (filtered_img + epsilon)).^2) ./ (kappa * q0_squared + epsilon));
        c_west = exp(-((west ./ (filtered_img + epsilon)).^2) ./ (kappa * q0_squared + epsilon));
        c_east = exp(-((east ./ (filtered_img + epsilon)).^2) ./ (kappa * q0_squared + epsilon));

        diffusion = lambda * (c_north .* north + c_south .* south + c_west .* west + c_east .* east);
        filtered_img = filtered_img + diffusion;
    end
end
