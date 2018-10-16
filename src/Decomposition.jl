module Decomposition

include("thirdparty/FastConv/FastConv.jl")
using DSP

export dwt

function morlet(F::T, fs::T) where {T}
    wtime = range(-1, stop=1, length=Int(2*fs))
    s = 6 / (2 * pi * F)
    wavelet = exp.(2 * 1im * pi * F * wtime) .* exp.(-wtime.^2 / (2*s^2))
    return wavelet |> Array{Complex{T}}
end

function dwt(data::Array{D, 2}, fs::T, frange::Array{T,1}; reflection::Bool=false, wavelet=morlet) where {D, T}
    data_fft = data

    if reflection
        data_flip = reverse(data_fft, dims=2)
        data_fft = [data_flip data_fft data_flip]
    end

    data_fft = data_fft |> Array{Complex{T}}
    result = zeros(Complex{T}, length(frange), size(data, 1), size(data, 2))

    for (fidx, ftarget) in enumerate(frange)
        w = wavelet(ftarget, fs) |> (x)->reshape(x, (1, length(x)))
        if reflection
            result[fidx, :, :] = conv2(data_fft, w)[:, Int(fs)+size(data, 2):end-Int(fs)-size(data, 2)]
        else
            result[fidx, :, :] = conv2(data_fft, w)[:, Int(fs):end-Int(fs)]
        end
    end
    return result

end

function dwt(data::Array{D, 1}, fs::T, frange::Array{T,1}; reflection::Bool=false, wavelet=morlet) where {D, T}
    data_fft = data

    if reflection
        data_flip = reverse(data_fft, dims=1)
        data_fft = [data_flip;data_fft;data_flip]
    end

    data_fft = data_fft |> Array{Complex{T}}
    result = zeros(Complex{T}, length(frange), size(data, 1))

    for (fidx, ftarget) in enumerate(frange)
        w = wavelet(ftarget, fs)
        if reflection
            result[fidx, :] = DSP.conv(data_fft, w)[Int(fs)+size(data, 1):end-Int(fs)-size(data, 1)]
        else
            result[fidx, :] = DSP.conv(data_fft, w)[Int(fs):end-Int(fs)]
        end
    end
    return result

end

end
