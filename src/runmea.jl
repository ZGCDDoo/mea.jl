push!(LOAD_PATH, "/home/hebertc1/Installations/Julia")
using Mea
using JSON

paramsfile = "statsparams0.json"
beta_ = JSON.parsefile(paramsfile)["beta"][1]

maxevals = 300_000
cutoff = 100.0

for ii in 0:10

try
	fname = "self_ctow" * string(ii) * ".dat"
	if !isfile(fname)
		break
	end
	foutdos = "dosjulia" * string(ii) * ".dat"
	foutdoscuba = "dosjulia_cuba" * string(ii) * ".dat"
	foutdostrace = "dosjulia_trace" * string(ii) * ".dat"
	foutdoscum = "dosjulia_cum" * string(ii) * ".dat"

	modelvec = Mea.Periodize.buildmodelvec(fname, paramsfile)

	coefs = Mea.Transport.coefstrans(modelvec, beta_, cutoff=cutoff, fout_name=foutdos, maxevals=maxevals)
	coefs_cuba = Mea.Transport.coefstrans(modelvec, beta_, cutoff=cutoff, fout_name=foutdoscuba, maxevals=maxevals, libintegrator="cuba")
	coefs_trace = Mea.Transport.coefstrans(modelvec, beta_, cutoff=cutoff, fout_name=foutdostrace, maxevals=maxevals,
										   libintegrator="cubature", fctper="make_akw2trace")
	coefs_cum =  Mea.Transport.coefstrans(modelvec, beta_, cutoff=cutoff, fout_name=foutdoscum, maxevals=maxevals,
											libintegrator="cubature", fctper="make_akw2cum")

	println("\n iteration = ", ii)
	println("coefs = ", coefs)
	println("coefs_cuba = ", coefs_cuba)
	println("coefs_trace = ", coefs_trace)
	println("coefs_cum = ", coefs_cum)

catch
end

end
