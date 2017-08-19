push!(LOAD_PATH, "/home/hebertc1/Installations/Julia")
using Mea
using JSON

paramsfile = "statsparams0.json"
beta = JSON.parsefile(paramsfile)["beta"][1]

for ii in 0:10

try
	fname = "self_ctow" * string(ii) * ".dat"
	if !isfile(fname)
		break
	end
	foutdos = "dosjulia" * string(ii) * ".dat"
	modelvec = Mea.Periodize.buildmodelvec(fname, paramsfile)
	coefs = Mea.Transport.coefstrans(modelvec, beta, cutoff=20.0, fout_name=foutdos, maxevals=200000)
	println("\n iteration = ", ii, "\ncoefstrans = \n")
	println(coefs)
catch
end

end
