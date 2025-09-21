function output = guided_filter(img, radius, epsilon)
    % Apply Guided Filter
    % img      : Input noisy image
    % radius   : Neighborhood size (radius of the filter)
    % epsilon  : Regularization parameter to smoothen output

    % Smooth the input image using the guided filter
    output = imguidedfilter(img, ...
                            'NeighborhoodSize', [radius radius], ...
                            'DegreeOfSmoothing', epsilon);
end
