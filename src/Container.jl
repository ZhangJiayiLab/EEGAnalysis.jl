


function create_epoch_bymarker(data::Array{D, 2}, marker::Array{T, 1}, 
        roi::Tuple{T,T}, fs::T, mbias::T=0.0) where {D, T}
    gap = (roi[2] - roi[1]) * fs |> ceil |> Int
    result = zeros(D, size(data, 1), gap, length(marker))
    for (midx, eachm) in enumerate(marker)
        start = (eachm + roi[1] + mbias) * fs |> floor |> Int
        result[:, :, midx] = data[:, start:start+gap-1]
    end
    return result
end

function create_epoch_bymarker(data::Array{D, 1}, marker::Array{T, 1}, 
        roi::Tuple{T,T}, fs::T, mbias::T=0.0) where {D, T}
    gap = (roi[2] - roi[1]) * fs |> ceil |> Int
    result = zeros(D, length(marker), gap)
    for (midx, eachm) in enumerate(marker)
        start = (eachm + roi[1] + mbias) * fs |> floor |> Int
        result[midx,:] = data[start:start+gap-1]
    end
    return result
end


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



########## ########## ########## ########## ##########
########## iSplit Data Container ########## ##########
########## ########## ########## ########## ##########

struct iSplitUnit
    chidx::Int
    edfname::String
    value::Array{Int16, 1}
    physicalunit::Float32
    samplingfrequency::Float32
end

struct iSplitContainer
    chidx::Int
    data::Dict{String, iSplitUnit}
end

function loadisplit(datadir, chname)
    chmatname = joinpath(datadir, @sprintf("sgch_ch%03d.mat", chname))
    raw = matread(chmatname)
    
    data = [name => iSplitUnit(chname, name, raw["values"][idx][:], 
                                raw["physicalunit"][idx], raw["samplingfrequency"][idx]) 
        for (idx,name) in enumerate(raw["edfnames"])] |> Dict
    
    return iSplitContainer(chname, data)
end


function chunk_isplit(isplitdata::iSplitContainer, markers::DataFrame, markername::String,
        roi::Tuple{T,T}, mbias::T=0; merge::Bool=false) where {T<:Float64}
    
    _result = Array{Array{Int16, 2}, 1}()
    
    for eachname in keys(isplitdata.data)
       markerlist = markers[(markers.id .== eachname).&(markers.mname .== markername), 3]
        
        _temp = create_epoch_bymarker(isplitdata.data[eachname].value, Array{T}(markerlist), 
                                      roi, T(isplitdata.data[eachname].samplingfrequency), mbias)
        
        push!(_result, _temp)
    end
    
    if merge
        result = _result[1]
        for i=2:length(_result)
            result = vcat(result, _result[i])
        end
    else
        result = _result
    end
    return result
end