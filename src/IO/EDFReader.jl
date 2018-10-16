module EDFReader

export readEDFFile, EDFData

struct EDFData
    
    # Header
    version::Int64                      # 8 ascii : version of this data format (0) 
    patient_info::String                # 80 ascii : local patient identification (mind item 3 of the additional EDF+ specs)
    record_info::String                 # 80 ascii : local recording identification (mind item 4 of the additional EDF+ specs)
    start_date::String                  # 8 ascii : startdate of recording (dd.mm.yy) (mind item 2 of the additional EDF+ specs)
    start_time::String                  # 8 ascii : starttime of recording (hh.mm.ss) 
    header_length::Int16                # 8 ascii : number of bytes in header record 
                                        # 44 ascii : reserved 
    recordnum::Int64                    # 8 ascii : number of data records (-1 if unknown, obey item 10 of the additional EDF+ specs) 
    sampleduration::Float32             # 8 ascii : duration of a data record, in seconds 
    nchannel::Int64                     # 4 ascii : number of signals (ns) in data record 
    channelLabels::Array{String, 1}     # ns * 16 ascii : ns * label (e.g. EEG Fpz-Cz or Body temp) (mind item 9 of the additional EDF+ specs)
    channelType::Array{String, 1}       # ns * 80 ascii : ns * transducer type (e.g. AgAgCl electrode) 
    physical_dim::Array{String, 1}      # ns * 8 ascii : ns * physical dimension (e.g. uV or degreeC) 
    physical_min::Array{Float32, 1}     # ns * 8 ascii : ns * physical minimum (e.g. -500 or 34) 
    physical_max::Array{Float32, 1}     # ns * 8 ascii : ns * physical maximum (e.g. 500 or 40) 
    digital_min::Array{Int16, 1}        # ns * 8 ascii : ns * digital minimum (e.g. -2048) 
    digital_max::Array{Int16, 1}        # ns * 8 ascii : ns * digital maximum (e.g. 2047) 
    prefiltering::Array{String, 1}      # ns * 80 ascii : ns * prefiltering (e.g. HP:0.1Hz LP:75Hz) 
    samples::Array{Int16, 1}            # ns * 8 ascii : ns * nr of samples in each data record 
    reserved_samples::Array{Int16, 1}   # ns * 32 ascii : ns * reserved
    
    # Data
    data::Array{Int16, 2}
    reserved_data::Array{Int16, 2}
    
    # meta
    Fs::Float32  #sampling frequency
    tspec::Array{Float32, 1}
    physical_unit::Array{Float32, 1}
end

function readEDFFile(filename::String)
    rawfile = open(filename)
    #Header
    version =      read!(rawfile, Array{UInt8, 1}(undef, 8))  .|> Char |> String |> Meta.parse
    patient_info = read!(rawfile, Array{UInt8, 1}(undef, 80)) .|> Char |> String |> strip
    record_info  = read!(rawfile, Array{UInt8, 1}(undef, 80)) .|> Char |> String |> strip
    start_date =   read!(rawfile, Array{UInt8, 1}(undef, 8))  .|> Char |> String |> strip
    start_time =   read!(rawfile, Array{UInt8, 1}(undef, 8))  .|> Char |> String |> strip
    header_length = read!(rawfile, Array{UInt8, 1}(undef, 8)) .|> Char |> String |> Meta.parse
    reserved = read!(rawfile, Array{UInt8, 1}(undef, 44))
    recordnum = read!(rawfile, Array{UInt8, 1}(undef, 8))     .|> Char |> String |> Meta.parse
    sampleduration = read!(rawfile, Array{UInt8, 1}(undef, 8)).|> Char |> String |> Meta.parse
    nchannel = read!(rawfile, Array{UInt8, 1}(undef, 4))      .|> Char |> String |> Meta.parse
    
    channelLabels = Array{String,1}(undef, nchannel)
    for chidx = 1:nchannel
        channelLabels[chidx] = read!(rawfile, Array{UInt8, 1}(undef, 16)) .|> Char |> String |> strip
    end

    channelType = Array{String,1}(undef, nchannel)
    for chidx = 1:nchannel
        channelType[chidx] = read!(rawfile, Array{UInt8, 1}(undef, 80)) .|> Char |> String |> strip
    end
    
    physical_dim = Array{String,1}(undef, nchannel)
    for chidx = 1:nchannel
        physical_dim[chidx] = read!(rawfile, Array{UInt8, 1}(undef, 8)) .|> Char |> String |> strip
    end
    
    physical_min = zeros(Float32, nchannel)
    for chidx = 1:nchannel
        physical_min[chidx] = read!(rawfile, Array{UInt8, 1}(undef, 8)) .|> Char |> String |> Meta.parse
    end
    
    physical_max = zeros(Float32, nchannel)
    for chidx = 1:nchannel
        physical_max[chidx] = read!(rawfile, Array{UInt8, 1}(undef, 8)) .|> Char |> String |> Meta.parse
    end
    
    digital_min = zeros(Int16, nchannel)
    for chidx = 1:nchannel
        digital_min[chidx] = read!(rawfile, Array{UInt8, 1}(undef, 8)) .|> Char |> String |> Meta.parse
    end
    
    digital_max = zeros(Int16, nchannel)
    for chidx = 1:nchannel
        digital_max[chidx] = read!(rawfile, Array{UInt8, 1}(undef, 8)) .|> Char |> String |> Meta.parse
    end
    
    prefiltering = Array{String, 1}(undef, nchannel)
    for chidx = 1:nchannel
        prefiltering[chidx] = read!(rawfile, Array{UInt8, 1}(undef, 80)) .|> Char |> String |> strip
    end
    
    samples = Array{Int16, 1}(undef, nchannel)
    for chidx = 1:nchannel
        samples[chidx] = read!(rawfile, Array{UInt8, 1}(undef, 8)) .|> Char |> String |> Meta.parse
    end
    
    reserved_samples = Array{Int16, 1}(undef, nchannel)
    for chidx = 1:nchannel
        reserved_samples[chidx] = read!(rawfile, Array{UInt8, 1}(undef, 32)) .|> Char |> String |> Meta.parse |> (x)-> (x==nothing) ? 0 : x
    end
    
    
    #Data
    data = zeros(Int16, (nchannel, recordnum*samples[1]))
    reserved_data = zeros(Int16, (nchannel, recordnum*reserved_samples[1]))
    
    step = samples |> sum
    for ri = 1:recordnum
        record_data = read!(rawfile, Array{Int16, 1}(undef, step))
        data[:, samples[1]*(ri-1)+1:samples[1]*ri] = reshape(record_data, (nchannel, samples[1]))
    end
    
    close(rawfile)
    
    Fs = samples[1] / sampleduration
    tspec = range(0, length=recordnum*samples[1], step=1/Fs) |> Array{Float32, 1}
    physical_unit = (physical_max .- physical_min) ./ (digital_max - digital_min)

    
    return EDFData( version, patient_info, record_info,
                    start_date, start_time, header_length,
                    recordnum, sampleduration, nchannel,
                    channelLabels, channelType, physical_dim, 
                    physical_min, physical_max, digital_min, digital_max,
                    prefiltering, samples, reserved_samples,
                    data, reserved_data, Fs, tspec, physical_unit)
end

end  # module EDFReader
