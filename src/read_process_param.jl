mutable struct proc_params

    rank::Int64
    niter::Int32
    f1::Float32
    f2::Float32

    c1::Float32
    alpha::Float32
    rho::Float32

    nt::Int32
    dt::Float32
    nx1::Int32
    nx2::Int32
end
 


function read_process_param(filename)

f = open(filename,"r")
fstring = read(f,String);
close(f)


ini = first(something(findlast("rank=",fstring),0:-1))
rank = parse(Float32, fstring[something(findnext(r"rank=.*",fstring,ini))][6:end]) 

ini = first(something(findlast("niter=",fstring),0:-1))
niter = parse(Float32, fstring[something(findnext(r"niter=.*",fstring,ini))][7:end])

ini = first(something(findlast("f1=",fstring),0:-1))
f1 = parse(Float32, fstring[something(findnext(r"f1=.*",fstring,ini))][4:end]) 

ini = first(something(findlast("f2=",fstring),0:-1))
f2 = parse(Float32, fstring[something(findnext(r"f2=.*",fstring,ini))][4:end]) 

ini = first(something(findlast("c1=",fstring),0:-1))
c1 = parse(Float32, fstring[something(findnext(r"c1=.*",fstring,ini))][4:end])

ini = first(something(findlast("alpha=",fstring),0:-1))
alpha = parse(Float32, fstring[something(findnext(r"alpha=.*",fstring,ini))][7:end])

ini = first(something(findlast("rho=",fstring),0:-1))
rho = parse(Float32, fstring[something(findnext(r"rho=.*",fstring,ini))][5:end])

nt=0
dt=0
nx1=0
nx2=0

params = proc_params(rank,niter,f1,f2,c1,alpha,rho,nt,dt,nx1,nx2)

return params

end