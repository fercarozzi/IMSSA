function doit(filein,fileout,pm,dm)

#>>>>>>>>>>>>>>>>> Read data  <<<<<<<<<<<<<<<<<<<<<

 di,hi,ext=SeisRead(filein);

 ni=size(di,2)
 xo=zeros(ni,2)


 xo[:,1] = SeisMain.ExtractHeader(hi,"sx")
 xo[:,2] = SeisMain.ExtractHeader(hi,"sy")

#>>>>>>>>>>>Define regular grid <<<<<<<<<<<<<<<<<<<<<

 a1_o = dm.minx 
 a1_f = dm.maxx

 a2_o = dm.miny
 a2_f = dm.maxy

 a1 = collect(a1_o:dm.dx:a1_f);

 a2 = collect(a2_o:dm.dy:a2_f)

 pm.nx1=size(a1,1)
 pm.nx2=size(a2,1)
 pm.dt=ext.d1
 pm.nt=ext.n1

#>>>>>>>> interpolation with regularization <<<<<<<<<<<<<<<<<<<

d_est = imssa(di,hi,ext;pm=pm,a1=a1,a2=a2,xo=xo,dx=dm.dx,dy=dm.dy);

xr=zeros(pm.nx1*pm.nx2,2)
aux=1

for j=1:pm.nx1
    for i=1:pm.nx2
    	xr[aux,1]=a1[j]
    	xr[aux,2]=a2[i]
	aux +=1
    end
end


#>>>>>>>> write to file <<<<<<<<<<<<<<<<<<<
writetofile(fileout,d_est,pm,dm,xr)

end