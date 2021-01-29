function test()

filein = "data/xs"

fileout = "data/imssa_xs"

data_param_file="grid_output.txt"

data_proc_file = "synth_param_proc.txt"

data_param = read_data_param(data_param_file)

proc_param = read_process_param(data_proc_file)

return filein,fileout,proc_param,data_param

end


