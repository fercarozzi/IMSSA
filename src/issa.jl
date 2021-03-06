"""
    issa(fi,fo,xi,xr1,xr2,dx1,dx2)

Reconstruction of irregular seismic data. 


# Output
* `f`: evaluated cost function at observed coordinates. 
* `g`: gradient at regular coordinates
"""
function issa(Dn,D0,xi,a1,a2,dx1,dx2)

#Dn is the result from the previous iteration in the w,x,y domain. It is regular
#D0 is the irregular, field data in the w,x,y domain.



# >>>>>>>>> Start reconstruction with f = S'*fi (adjoint solution) <<<<<<<<<<<<<

 
 fi = zero(D0);
 
# apply interpolator to go to the irregular coordinates.

for i=1:size(Dn,1)
    fi[i,:] = interpol2d(Dn[i,:,:],false;xi=xi,xr1=a1,xr2=a2,dx1=dx1,dx2=dx2)
end



# calculate the error between the estimation and the initial data, in irregular coord
       f = fi-D0;

#bring the error to regular coords. Allows calculation of next iteration x_{n+1}

g = zero(Dn);
for i =1 :size(Dn,1)
    g[i,:,:] = interpol2d(f[i,:],true;xi=xi,xr1=a1,xr2=a2,dx1=dx1,dx2=dx2)
end

return norm(f),g
end


#***************************************************************

"""
    fx(Dn,D0,xi,a1,a2,dx1,dx2)

Evaluates cost function at observed coordinates. 
"""
function fx(Dn,D0,xi,a1,a2,dx1,dx2)

#Dn is the result from the previous iteration in the w,x,y domain. It is regular
#D0 is the irregular, field data in the w,x,y domain.

 
 fi = zero(D0); 
 
# apply interpolator to go to the irregular coordinates.

for i=1:size(Dn,1)
    fi[i,:] = interpol2d(Dn[i,:,:],false;xi=xi,xr1=a1,xr2=a2,dx1=dx1,dx2=dx2)
end

# calculate the error between the estimation and the initial data, in irregular coord
       f = fi-D0;

return norm(f)
end



