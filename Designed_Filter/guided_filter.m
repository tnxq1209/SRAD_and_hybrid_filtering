function guided_img = guided_filter(I, P, radius, epsilon)
    % I: Guidance image (input image)
    % P: Filtering image (input image)
    % radius: Window radius
    % epsilon: Regularization parameter

    [rows, cols] = size(I);
    N = box_filter(ones(rows, cols), radius);

    mean_I = box_filter(I, radius) ./ N;
    mean_P = box_filter(P, radius) ./ N;
    corr_I = box_filter(I .* I, radius) ./ N;
    corr_IP = box_filter(I .* P, radius) ./ N;

    var_I = corr_I - mean_I .* mean_I;
    cov_IP = corr_IP - mean_I .* mean_P;

    a = cov_IP ./ (var_I + epsilon);
    b = mean_P - a .* mean_I;

    mean_a = box_filter(a, radius) ./ N;
    mean_b = box_filter(b, radius) ./ N;

    guided_img = mean_a .* I + mean_b;
end

function result = box_filter(img, radius)
    % Fast box filter for mean calculation
    kernel = ones(2*radius + 1, 2*radius + 1);
    result = conv2(img, kernel, 'same');
end
