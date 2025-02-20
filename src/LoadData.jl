include(scriptsdir("Params.jl"))

ffrf(file) = fill_from_root_file(file, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"])

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

function load_data_processes(dir::String, mode::String)
    filesDict = Dict()
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
            push!(
                processes,
                DataProcess(
                    df.reconstructedEnergy1 .+ df.reconstructedEnergy2,
                    sumEParams[Symbol(fileName)]
                )
            )
        elseif mode == "singleE"
            push!(
                processes,
                DataProcess(
                    vcat(df.reconstructedEnergy1, df.reconstructedEnergy2),
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
    filesDict = Dict()
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




function load_ndim_processes(dir::String, binsAngle, binsESingle, binsESum)
    filesDict = Dict()
    full_dir = datadir("sims", dir)
    processes = NDimDataProcess[]

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

        df = ffrf(f)
        fileName = split(file, ".")[1]  |> split |> first 
        
        push!(
                processes,
                NDimDataProcess(
                    df.phi,
                    [maximum([df.reconstructedEnergy1[i], df.reconstructedEnergy2[i]]) for i in 1:nrow(df)],
                    [sum([df.reconstructedEnergy1[i], df.reconstructedEnergy2[i]]) for i in 1:nrow(df)],
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
