using Mea
using Base.Test

tests = [
        "green",
        #"periodize",
        "sigmadc"
        ]


if length(ARGS) > 0
    tests = ARGS
end


for t in tests
    test_file = "test_$t.jl"
    print_with_color(:green, "* $test_file\n")
    include(test_file)
end
