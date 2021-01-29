"""
	writetofile(fout,d,pm,dm,xr)
Writes output to file with SeismicJulia format

# Arguments
* `fout` : Output file name
* `d` : Reconstructed traces
* `pm`: processing parameters
* `dm`: data parameters
* `xr`: Regular coordinates
"""
function writetofile(fileout,d,pm,dm,xr) 

  ntrace=size(d,3)*size(d,2)

  h = Array{Header}(undef,ntrace);
  j=1


  for ix1 = 1:size(d,2)
    for ix2 = 1:size(d,3)
      h[j] = SeisMain.InitSeisHeader();
      h[j].tracenum = j;
      h[j].d1 = pm.dt;
      h[j].n1 = pm.nt;
      h[j].sy = dm.minx+dm.dx*(ix1-1);
      h[j].sx = dm.miny+dm.dy*(ix2-1);
      
      j +=1
     end
  end


 extent = SeisMain.Extent(convert(Int32, size(d,1)), convert(Int32, size(d,2)),
                            convert(Int32, size(d,3)), 1,1,
			    convert(Float32, 0),
                            convert(Float32, minimum(xr[:,2])),
			    convert(Float32, minimum(xr[:,1])),
                            convert(Float32, 1),
			    convert(Float32, 1),
                            convert(Float32, pm.dt), convert(Float32, dm.dx),
                            convert(Float32, dm.dy), convert(Float32, 1),
                            convert(Float32,1), "Time", "sx", "sy", "gx",
                            "gy", "s", "m", "m", "m", "m", "")


SeisWrite(fileout,d,h,extent)

end
