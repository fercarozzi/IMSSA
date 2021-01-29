using SeisMain,FFTW,DSP,LinearAlgebra,TSVD,SeisProcessing
using PyPlot, SeisPlot
using Statistics

 include("src/read_data_param.jl")
 include("src/read_process_param.jl")

 include("src/backtracking.jl")
 include("src/doit.jl")
 include("src/FFTfunc.jl")
 include("src/imssa.jl")
 include("src/issa.jl")
 include("src/operator.jl")
 include("src/ssa3d.jl")
 include("src/writetofile.jl")
 
 include("tests.jl")

 filein,fileout,proc_param,data_param = test()

 println("input file=", filein," output file= ",fileout)


 doit(filein,fileout,proc_param,data_param)


d0,h0,e0 = SeisRead("data/reg_xs")
dout,hout,eout = SeisRead("data/imssa_xs")

Q = 10*log(10,sum(d0.^2)/sum((d0-dout).^2))

println("S/N_o reconstrution= ",Q)

#--------------------------------------------
# Plot
figure(1,(10,10))

subplot(121)
SeisPlotTX(d0[:,5,:],style="wiggles",xlabel="X coordinate (m)",ox=1,dx=data_param.dx,dy=e0.d1,ylabel="Time (s)",fignum=1,title="Ideal")

subplot(122)
SeisPlotTX(dout[:,5,:],style="wiggles",xlabel="X coordinate (m)",ox=1,dx=data_param.dx,dy=e0.d1,ylabel="Time(s)",fignum=1,title="Reconstruction")
tight_layout()
