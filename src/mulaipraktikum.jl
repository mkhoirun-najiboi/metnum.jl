function mulaipraktikum()
    repo_directory = joinpath(@__DIR__,".julia","packages","metnum")
    folders = readdir(repo_directory)
    unix_time = Array{Any}(undef,length(folders))
    k=1
    for folder in folders
        date_time = unix2datetime(mtime(joinpath(repo_directory, folder)))
        unix_time[k] = Int(floor(datetime2unix(date_time)))
        k+=1
    end
    path = joinpath(repo_directory, folders[argmax(unix_time)],"notebookpraktikum")
    IJulia.notebook(;dir=path)
end
