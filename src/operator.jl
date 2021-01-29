"""
    interpol2d (in,adj;keywords)

Interpolator operator for 3D reconstruction of irregular seismic data

# Arguments
* `in` : input data
* `adj::Bool`: true for operator, false for adjoint

# keywords
* `xi::Array{Float64,2}`: matrix with irregular coordinates
* `xr1,xr2::Array{Float64,1},Array{Float64,1}`: 2 vectors with regular coordinates
* `dx1`,`dx2` : sampling rate

# Output
* `out`: output data
"""
function interpol2d(in,adj;xi=1,xr1=1,xr2=1,dx1=1,dx2=1)
ni = size(xi,1)
nr1= size(xr1,1)
nr2= size(xr2,1)

if (adj)

    out=complex(zeros(nr1,nr2))


    for k=1:ni
        ia = convert(Int64,floor(abs(xi[k,1]-xr1[1])/dx1)+1)
	ja = convert(Int64,floor(abs(xi[k,2]-xr2[1])/dx2)+1)

	ib = ia + 1
	jb = ja + 1


	if  (1<ib<=nr1 && 1<jb<=nr2)
		t = abs((xi[k,1]-xr1[ia]))/dx1
		u = abs((xi[k,2]-xr2[ja]))/dx2

		out[ia,ja] += (1-t) * (1-u) * in[k]
		out[ib,ja] +=    t  * (1-u) * in[k]
		out[ib,jb] +=    t  *    u  * in[k]
		out[ia,jb] += (1-t) *    u  * in[k]
         end
     end
 else
    out=complex(zeros(ni))
    for k=1:ni
        ia = convert(Int64,floor(abs(xi[k,1]-xr1[1])/dx1)+1)
    	ja = convert(Int64,floor(abs(xi[k,2]-xr2[1])/dx2)+1)

    	ib = ia + 1
    	jb = ja + 1

    	if  (1<ib<=nr1 && 1<jb<=nr2)
    	   	t = abs((xi[k,1]-xr1[ia]))/dx1
    		u = abs((xi[k,2]-xr2[ja]))/dx2

	        out[k] += (1-t) * (1-u) * in[ia,ja] +
		  	             t  * (1-u) * in[ib,ja] +
			             t  *    u  * in[ib,jb] +
			          (1-t) *    u  * in[ia,jb]

        end
    end

 end

return out

end
