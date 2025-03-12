# include(scriptsdir("Params.jl"))
using Distributions

ffrf(file) = fill_from_root_file(file, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"])
function smear_energy(E::Real, fwhm::Real) 
    sigma = fwhm/2.355 * E * sqrt(1000/E)
    return rand(Normal(E, sigma))
end

function load_files(dir::String)
    filesDict = Dict()
    full_dir = datadir("sims", dir)

    println("Loading files from: $full_dir ...")

    nFiles = 0
    for file in readdir(full_dir)
        nFiles += 1
        df = ROOTFile(joinpath(full_dir, file)) |> ffrf
        fileName = split(file, ".")[1]  |> split |> first 

        filesDict[fileName] = df
    end
    println("Loaded $nFiles files.")

    return filesDict
end

function load_data_processes(dir::String, mode::String; fwhm = 0.08)
    include(scriptsdir("Params.jl"))
    full_dir = datadir("sims", dir)
    processes = DataProcess[]

    println("Loading files from: $full_dir ...")
    println("mode: $mode ")

    nFiles = 0
    for file in readdir(full_dir)
        nFiles += 1
        if( split(file, ".")[end] != "root" )
            continue
        end
        
        f = ROOTFile(joinpath(full_dir, file)) 
        if(!haskey(f, "tree"))
            continue
        end

        df = ffrf(f)
        fileName = split(file, ".")[1]  |> split |> first 
        
        if mode == "sumE"
            e = df.reconstructedEnergy1 .+ df.reconstructedEnergy2
            push!(
                processes,
                DataProcess(
                    ifelse.(fwhm > 0, smear_energy.(e, fwhm), e) ,
                    sumEParams[Symbol(fileName)]
                )
            )
        elseif mode == "singleE"
            e = vcat(df.reconstructedEnergy1, df.reconstructedEnergy2)
            push!(
                processes,
                DataProcess(
                    ifelse.(fwhm > 0, smear_energy.(e, fwhm), e) ,
                    singleEParams[Symbol(fileName)]
                )
            )
        elseif mode == "phi"
            push!(
                processes,
                DataProcess(
                    df.phi,
                    phiParams[Symbol(fileName)]
                )
            )
        else
            @error "mode must be one of: sumE, singleE, phi. Chosen mode: $mode!"
        end
    end
    println("Loaded $nFiles files.")

    return processes
end



function load_hist_processes(dir::String, mode::String)
    include(scriptsdir("Params.jl"))
    full_dir = datadir("sims", dir)
    processes = HistProcess[]

    println("Loading files from: $full_dir ...")
    println("mode: $mode ")

    nFiles = 0
    for file in readdir(full_dir)
        nFiles += 1
        if( split(file, ".")[end] != "root" )
            continue
        end
        
        f = ROOTFile(joinpath(full_dir, file)) 
        if(haskey(f, "tree"))
            continue
        end
        fileName = split(file, ".")[1]  |> split |> first 

        for k in keys(f) # if there are multiple histograms in the root file we save each?
            h1 = UnROOT.parseTH(f[k], raw= false)
            if mode == "sumE"
                push!(
                    processes,
                    HistProcess(
                        h1,
                        k,
                        sumEParams[Symbol(fileName)]
                    )
                )

            elseif mode == "singleE"
                push!(
                    processes,
                    HistProcess(
                        h1,
                        k,
                        singleEParams[Symbol(fileName)]
                    )
                )
            elseif mode == "phi"
                push!(
                    processes,
                    HistProcess(
                        h1,
                        k,
                        phiParams[Symbol(fileName)]
                    )
                )
            else
                @error "mode must be one of: sumE, singleE, phi. Chosen mode: $mode!"
            end
        end
    end
    println("Loaded $nFiles files.")

    return processes
end




function load_ndim_processes(dir::String, bins::NamedTuple, varNames::Vector{String}; fwhm = 0.08)
    include(scriptsdir("Params.jl"))
    full_dir = datadir("sims", dir)
    processes = DataProcess3D[]

    println("Loading files from: $full_dir ...")
    println("mode: NDim ")

    nFiles = 0
    for file in readdir(full_dir)
        nFiles += 1
        if( split(file, ".")[end] != "root" )
            continue
        end
        
        f = ROOTFile(joinpath(full_dir, file)) 
        if(!haskey(f, "tree"))
            continue
        end

        data = LazyTree(f, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"])

        Ei = Vector{Float64}(undef, length(data)) # Float32[] # individual energies
        Es = Vector{Float64}(undef, length(data)) #Float32[] # sum of energies
        phi = Vector{Float64}(undef, length(data)) # Float32[] # phi angle

        Threads.@threads for (i,event) in enumerate(data)
            i%1_000_000 == 0 && println("$i/$(length(data))  events processed!")
            Ei_ = maximum([event.reconstructedEnergy1, event.reconstructedEnergy2])
            Es_ = sum([event.reconstructedEnergy1, event.reconstructedEnergy2])

            Ei[i] = smear_energy(Ei_, fwhm)
            phi[i] = event.phi
            Es[i] = smear_energy(Es_, fwhm)
        end

        fileName = split(file, ".")[1]  |> split |> first 
        
        push!(
                processes,
                DataProcess3D(
                    phi,
                    Ei,
                    Es,
                    binsAngle,
                    binsESingle,
                    binsESum,
                    singleEParams[Symbol(fileName)]
                )
            )
        println("$fileName loaded")

    end
    println("Loaded $nFiles files.")

    return processes
end

function load_3D_processes(dir::String, binsAngle, binsESingle, binsESum; fwhm = 0.08)
    include(scriptsdir("Params.jl"))
    full_dir = datadir("sims", dir)
    processes = DataProcess3D[]

    println("Loading files from: $full_dir ...")
    println("mode: NDim ")

    nFiles = 0
    for file in readdir(full_dir)
        nFiles += 1
        if( split(file, ".")[end] != "root" )
            continue
        end
        
        f = ROOTFile(joinpath(full_dir, file)) 
        if(!haskey(f, "tree"))
            continue
        end

        data = LazyTree(f, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"])

        Ei = Vector{Float64}(undef, length(data)) # Float32[] # individual energies
        Es = Vector{Float64}(undef, length(data)) #Float32[] # sum of energies
        phi = Vector{Float64}(undef, length(data)) # Float32[] # phi angle

        Threads.@threads for (i,event) in enumerate(data)
            i%1_000_000 == 0 && println("$i/$(length(data))  events processed!")
            Ei_ = maximum([event.reconstructedEnergy1, event.reconstructedEnergy2])
            Es_ = sum([event.reconstructedEnergy1, event.reconstructedEnergy2])

            Ei[i] = smear_energy(Ei_, fwhm)
            phi[i] = event.phi
            Es[i] = smear_energy(Es_, fwhm)
        end

        fileName = split(file, ".")[1]  |> split |> first 
        
        push!(
                processes,
                DataProcess3D(
                    phi,
                    Ei,
                    Es,
                    binsAngle,
                    binsESingle,
                    binsESum,
                    singleEParams[Symbol(fileName)]
                )
            )
        println("$fileName loaded")

    end
    println("Loaded $nFiles files.")

    return processes
end