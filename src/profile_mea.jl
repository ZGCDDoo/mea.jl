push!(LOAD_PATH, pwd())
using Mea.Periodize
using Mea.Transport
using ProfileView


function code_to_profile()
modelvec = Periodize.buildmodelvec("self_ctow0.dat", "statsparams0.json")
Periodize.calcdos2(modelvec)
end

#taken from: https://thirld.com/blog/2015/05/30/julia-profiling-cheat-sheet/
function benchmark()
    # Any setup code goes here.

    # Run once, to force compilation.
    println("======================= First run:")
    srand(69)
    @time code_to_profile()

    # Run a second time, with profiling.
    println("\n\n======================= Second run:")
    srand(69)
    Profile.init(delay=0.01)
    Profile.clear()
    Profile.clear_malloc_data()
    @profile @time code_to_profile()

    # Write profile results to profile.bin.
    r = Profile.retrieve()
    f = open("profile.bin", "w")
    serialize(f, r)
    close(f)
end


function view_results()
using ProfileView
f=open("profile.bin")
r=deserialize(f);
ProfileView.view(r[1], lidict=r[2])

end