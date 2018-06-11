using PlasmoAlgorithms
using Plasmo

function genproblem(sublib::String,params...)
    numparams = length(params)
    basename = sublib
    for i in 1:numparams
        basename *= "_$(params[i])"
    end
    g = smpsread("$sublib/$basename")
    return g
end
