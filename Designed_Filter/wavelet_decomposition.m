function wavelet_img = wavelet_decomposition(img, wavelet_type, level)
    if level >= 1
        [cA1, cH1, cV1, cD1] = dwt2(img, wavelet_type);
    end

    if level >= 2
        [cA2, cH2, cV2, cD2] = dwt2(cA1, wavelet_type);
        cA1_reconstructed = idwt2(cA2, cH2, cV2, cD2, wavelet_type);
    else
        cA1_reconstructed = cA1;
    end

    wavelet_img = idwt2(cA1_reconstructed, cH1, cV1, cD1, wavelet_type);
end
