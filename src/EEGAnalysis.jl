module EEGAnalysis

using MAT

include("Container.jl")
include("Decomposition.jl")
include("IO/EDFReader.jl")

export Container
export Decomposition
export EDFReader
export EDFData

end  # module EEGAnalysis
