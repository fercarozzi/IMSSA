"""
    ssa_3d(d,p)
ssa_3d applies mssa filtering to a 3D cube.

# Arguments
* `d`: 2D array for one frequency slide.
* `p`: desired rank (number of dips)

# Output
* `dout`: filtered data

"""
function ssa_3d(d,p)

nx,ny = size(d);

# Size of the Hankel Matrix for y.
Ncol = convert(Int64,floor(ny/2))+1;
Nrow = ny-Ncol+1;

# Size of HAnkel MAtrix of Hankel Matrixes in x.
Lcol = convert(Int64,floor(nx/2))+1;
Lrow = nx-Lcol+1;


M = complex(zeros(Lrow*Nrow,Lcol*Ncol))
tmp_fx = complex(zeros(Lcol+Lrow-1,ny))

    for lc = 1:Lcol
        for lr = 1:Lrow

            tmp_fx[lr+lc-1,:]  = (d[lr+lc-1,:])';

            for ic = 1:Ncol;
                for ir = 1:Nrow;
		M[(lr*Nrow)-Nrow+ir,(lc*Ncol)-Ncol+ic] = tmp_fx[lr+lc-1,ir+ic-1];
                end
            end

        end
    end


# SVD decomposition with P largest sing values 

U,S,V = tsvd(M,p);
Mout = U*Diagonal(S)*V'

# Average along anti-diagonals

    dout=complex(zero(d));

    Count = zeros(ny,nx);
    tmp = complex(zeros(ny,nx));
    tmp2 = complex(zeros(ny,nx));

    for lc = 1:Lcol
        for lr = 1:Lrow
            for ic = 1:Ncol;
                for ir = 1:Nrow;
                    Count[ir+ic-1,lr+lc-1] += 1;
                    tmp[ir+ic-1,lr+lc-1]  += Mout[(lr*Nrow)-Nrow+ir,(lc*Ncol)-Ncol+ic];
                end;
            end


            tmp2[:,lr+lc-1] = tmp[:,lr+lc-1]./Count[:,lr+lc-1];

            dout[lr+lc-1,:] = tmp2[:,lr+lc-1]';

        end
    end



return dout
end
