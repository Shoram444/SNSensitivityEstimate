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

function load_processes(dir::String, mode::String)
    filesDict = Dict()
    full_dir = datadir("sims", dir)
    processes = []

    println("Loading files from: $full_dir ...")
    println("mode: $mode ")

    nFiles = 0
    for file in readdir(full_dir)
        nFiles += 1
        if( split(file, ".")[end] != "root" )
            continue
        end
        
        df = ROOTFile(joinpath(full_dir, file)) |> ffrf
        fileName = split(file, ".")[1]  |> split |> first 
        
        if mode == "sumE"
            push!(
                processes,
                Process(
                    df.reconstructedEnergy1 .+ df.reconstructedEnergy2,
                    sumEParams[Symbol(fileName)]
                )
            )
        elseif mode == "singleE"
            push!(
                processes,
                Process(
                    vcat(df.reconstructedEnergy1, df.reconstructedEnergy2),
                    singleEParams[Symbol(fileName)]
                )
            )
        elseif mode == "phi"
            push!(
                processes,
                Process(
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
