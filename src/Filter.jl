function gaussian_kernel(fs::T, sigma::T; dims=nothing) where {T}
    ktime = range(-1, stop=1, length=Int(2*fs))
    kernel = 1/(sqrt(2*pi)*sigma) * exp.(-ktime.^2 / (2*sigma^2))
    if dims == nothing
        return kernel |> Array{T}
    elseif dims == 1
        return kernel |> Array{T} |> x->reshape(x, (1, length(x)))
    elseif dims == 2
        return kernel |> Array{T} |> x->reshape(x, (length(x), 1))
    else
        return false
    end
end

function gaussianwind(data::Array{D, 1}, fs::T, sigma::T) where {D, T}
    k = gaussian_kernel(fs, sigma)
    totalpwr_filter = DSP.conv(data, k);
    return totalpwr_filter[Int(fs):end-Int(fs)] / sum(k)
end

function gaussianwind(data::Array{D, 2}, fs::T, sigma::T) where {D, T}
    k = gaussian_kernel(fs, sigma)
    totalpwr_filter = zeros(size(data))
    
    for idx = 1:size(data, 1)
        totalpwr_filter[idx, :] = gaussianwind(data[idx, :], fs, sigma)
    end
    totalpwr_filter
end
    