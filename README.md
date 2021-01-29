The package includes a demo of I-MSSA developed for julia1.1.

caveat: "If no such file or directory" error, generate the synthetic files again. 
To run execute 

> include("main.jl")


Dependencies to install before testing: 
- SeisMain: https://github.com/SeismicJulia/SeisMain.jl.git
- SeisProcessing: https://github.com/SeismicJulia/SeisProcessing.jl.git
- SeisPlot: https://github.com/SeismicJulia/SeisPlot.jl.git
 Registered packages:
- FFTW: 
- TSVD  
- LinearAlgebra
- PyPlot

folder /data has the input data
folder /data_process has the codes for data modeling
folder /src has the source codes for imssa.

If you use this package, please cite
Carozzi F. and Sacchi, M.D., 2021. Interpolated Multichannel Singular Spectrum Analysis (I-MSSA): a reconstruction method that honors true trace coordinates: Geophysics.