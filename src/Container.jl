module Container

import MAT

struct CompactDataContainer
    # meta
    name::String
    resultdir::String
    data::String

    # data
    channels::Array{Float64, 2}
    times::Array{Float64, 1}
    markers::Dict{String, Array{Float64, 1}}

    # constants
    fs::Int64
    roi::Tuple{Float64, Float64}
end

function createCompactDataContainer(datadir::String, resultdir::String,
    patientname::String, expname::String, fs::Int64, roi::Tuple{Float64, Float64})

    _resultdir = joinpath(resultdir, patientname)
    _datadir = joinpath(datadir, patientname, "EEG")

    #TODO: _checkdir(_resultdir, expname)

    rawdata = MAT.matread(joinpath(_datadir, "Compact", expname * ".mat"))

    channels = rawdata["channels"]
    times = rawdata["times"][:]
    markers = [i => [rawdata["markers"][i];][:] for i in keys(rawdata["markers"])] |> Dict

    # TODO: roi with iti


    return CompactDataContainer(expname, _resultdir, _datadir,
        channels, times, markers, fs, roi)
end

function createEpochByMarkername(container::CompactDataContainer, chidx::Int, markername::String)
    gap = Int((container.roi[2] - container.roi[1])*container.fs)
    marker = container.markers[markername]
    groupdata = zeros(Float64, length(marker), gap)
    for (midx, eachm) in enumerate(marker)
        starter = (eachm + container.roi[1]) * container.fs |> floor |> Int
        groupdata[midx, :] = container.channels[chidx, starter:starter+gap-1]
    end
    return groupdata
end

function createEpochByMarker(data::Array{D,1}, marker::Array{T,1}, roi::Tuple{T,T}, fs::T; mbias::T=0.0) where {T, D}
    gap = Int((roi[2] - roi[1])*fs)
    groupdata = zeros(D, length(marker), gap)
    for (midx, eachm) in enumerate(marker)
        starter = (eachm + mbias + roi[1]) * fs |> floor |> Int
        groupdata[midx, :] = data[starter:starter+gap-1]
    end
    return groupdata
end

function createEpochByMarker(data::Array{D,2}, marker::Array{T,1}, roi::Tuple{T,T}, fs::T; mbias::T=0.0) where {T, D}
    gap = Int((roi[2] - roi[1])*fs)
    groupdata = zeros(D, size(data, 1), length(marker), gap)
    for idx = 1:size(data, 1)
        groupdata[idx, :, :] = createEpochByMarker(data[idx, :], marker, roi, fs)
    end
    return groupdata
end

########## ########## ########## ########## ##########
########## Split Data Container  ########## ##########
########## ########## ########## ########## ##########

struct SplitDataContainer
    datadir::String
    Chidx::Int64

    fs::Int64
    #TODO: marker_bias::Dict{String, Float64}
    ERP::Dict{String, Array{Float64, 2}}
end



end
