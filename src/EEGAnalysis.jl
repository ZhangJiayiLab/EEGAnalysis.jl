module EEGAnalysis

using MAT
using DSP
using Printf
using DataFrames
using CSV

include("Container.jl")
include("Decomposition.jl")
include("IO/EDFReader.jl")
include("Filter.jl")

export loadcompact, loadedf, create_epoch_bymarker, chunk_isplit, loadisplit
export gaussianwind
export Decomposition


"""load European Data Format file

Syntax: EDFData = loadedf(filename)

Key argument:
filename -- (String) file path

Return:
EDFData  -- (EDFData)

Notes:

struct EDFData
    
    # Header
    version::Int64                      # version of this data format
    patient_info::String                # local patient identification
    record_info::String                 # local recording identification
    start_date::String                  # startdate of recording (dd.mm.yy)
    start_time::String                  # starttime of recording (hh.mm.ss) 
    header_length::Int16                # number of bytes in header record 
                                        # reserved 
    recordnum::Int64                    # number of data records
    sampleduration::Float32             # duration of a data record, in seconds 
    nchannel::Int64                     # number of signals (ns) in data record 
    channelLabels::Array{String, 1}     # ns * label (e.g. EEG Fpz-Cz or Body temp)
    channelType::Array{String, 1}       # ns * transducer type (e.g. AgAgCl electrode) 
    physical_dim::Array{String, 1}      # ns * physical dimension (e.g. uV or degreeC) 
    physical_min::Array{Float32, 1}     # ns * physical minimum (e.g. -500 or 34) 
    physical_max::Array{Float32, 1}     # ns * physical maximum (e.g. 500 or 40) 
    digital_min::Array{Int16, 1}        # ns * digital minimum (e.g. -2048) 
    digital_max::Array{Int16, 1}        # ns * digital maximum (e.g. 2047) 
    prefiltering::Array{String, 1}      # ns * prefiltering (e.g. HP:0.1Hz LP:75Hz) 
    samples::Array{Int16, 1}            # ns * nr of samples in each data record 
    reserved_samples::Array{Int16, 1}   # ns * reserved
    
    # Data
    data::Array{Int16, 2}
    reserved_data::Array{Int16, 2}
    
    # meta
    Fs::Float32  #sampling frequency
end

"""
function loadedf(filename::String)
    EDFReader.readEDFFile(filename)
end


"""

"""
function loadcompact(datadir::String, resultdir::String, patientname::String, 
        expname::String, fs::Int64, roi::Tuple{Float64, Float64})

    createCompactDataContainer(datadir::String, 
        resultdir::String, patientname::String, expname::String, 
        fs::Int64, roi::Tuple{Float64, Float64})
end
    
end  # module EEGAnalysis
