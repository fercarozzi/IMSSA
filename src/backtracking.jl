"""
	backtracking(d,dn,x,xr1,xr2,dx1,dx2; <keyword args>)
	
Find a step size for the steepest descent via line-search 

# Arguments
* `d`: 3D seismic data. 
* `dn`: Previous approximation.
* `x` : Irregular coordinates. 
* `xr1',`xr2` : Regular coordinates to form a mesh
* `dx1` `dx2`: Sampling interval

# Keyword Arguments
* Structure with processing parameters
"""
function BackTracking(Dn,D0,xo,a1,a2,dx,dy;pm=proc_params[])

# initialize values

     max_it = 500

     c_1 = pm.c1
     alpha = pm.alpha
     rho = pm.rho
    
    iteration = 0

    dx1 = dx
    dx2 = dy

    # Backtrack until it satisfies sufficient decrease condition

    f0,df0 = issa(Dn,D0,xo,a1,a2,dx1,dx2)

    x1 = Dn - alpha .* df0
    f1 = fx(x1,D0,xo,a1,a2,dx1,dx2)

    println("initial condition in backtrack: ",norm(f1)- norm(f0 + c_1 * alpha * dot(df0,-df0))," >0")

    while norm(f1)- norm(f0 + c_1 * alpha * dot(df0,-df0)) > 0

    	# Increment the number of steps 
        iteration += 1

        # Ensure termination
        if iteration > max_it
            println("Linesearch failed to converge, reached maximum iterations")
	    return
        end

        # Shrink proposed step-size:

	alpha = alpha * rho

        # Evaluate f(x) at proposed position

	x1 = Dn - alpha .* df0
    	f1 = fx(x1, D0,xo,a1,a2,dx1,dx2)

    end
 
    return alpha
end