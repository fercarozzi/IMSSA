"""
    imssa(d,h,ext;<keyword args>)

I-MSSA engine

# Input
* `d`,`h`,`e` : Field data with SeismicJulia format

# Keyword arguments
* `pm`::proc_params : Process parameters
* `a1`,`a2` : Regular coordinates to form a mesh
* `xo` : Observed coordinates
* `dx`,`dy` : Sampling interval
"""
function imssa(d,h,e;pm=proc_params[], a1=1,a2=1,xo=1,dx=1,dy=1)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Calculate FFT of observed traces

 D0, if1,if2,nf = FFTfunc(d,pm.dt,pm.f1,pm.f2);

 D0 = D0[if1:if2,:,:]
 
# Regularize observed traces

dr = zeros(pm.nt,pm.nx1,pm.nx2);

 for i =1:pm.nt
    dr[i,:,:] .= interpol2d(d[i,:],true; xi=xo,xr1=a1,xr2=a2,dx1=dx,dx2=dy)
 end

#dr = dr 

#  FFT of Regularized observed traces

 Dn,if1,if2,nf = FFTfunc(dr,pm.dt,pm.f1,pm.f2);

 Dn = Dn[if1:if2,:,:]
 
 nw = size(Dn,1)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Iterate to solve the objective function

Dp = zero(Dn);

df = zero(Dn);

for i =1:pm.niter
    println("  >> Iteration No ",i)

# Objective function step

   mu = BackTracking(Dn,D0,xo,a1,a2,dx,dy,pm=pm) 

   println("step size= ",mu) 

   df0 = zero(Dn);
   f0,df0 = issa(Dn,D0,xo,a1,a2,dx,dy);
   
# New approx 

   Dn1 = Dn - mu* df0;

# Project to subset of low rank matrices, D projected = Dp

 Dp = zero(Dn1);

 for k =1:nw
   Dp[k,:,:] = ssa_3d(Dn1[k,:,:],pm.rank)
 end

 Dn = Dp
end


#>>>>>>>>>>>>> Antitransform to time domain, calculate symmetries<<<

  D_est = complex(zeros(nf,pm.nx1,pm.nx2));
  
  D_est[if1:if2,:,:] .= Dp

for k =if1:if2
    D_est[nf-k+2,:,:]= conj(D_est[k,:,:])
end

 d_est = ifft(D_est,1);
 d_est = real(d_est[1:pm.nt,:,:]);

return d_est
end 
