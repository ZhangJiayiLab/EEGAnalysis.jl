function gaussian_kernel(fs::T, sigma::T) where {T}
    ktime = range(-1, stop=1, length=Int(2*fs))
    kernel = 1/(sqrt(2*pi)*sigma) * exp.(-ktime.^2 / (2*sigma^2))
    return kernel |> Array{T}
end

function gaussianwind(data::Array{D, 1}, fs::T, sigma::T) where {D, T}
    totalpwr_filter = DSP.conv(data, gaussian_kernel(fs, sigma));
    return totalpwr_filter[fs:end-fs]
end