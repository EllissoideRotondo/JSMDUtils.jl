module FileUtils

using JSON3: JSON3

import YAML as YAMLlib
using JSMDInterfaces.FilesIO

for fn in (:JSON, :TXT, :YAML)
    @eval begin
        FilesIO.@filetype $fn FilesIO.AbstractFile
        export $fn
    end
end

"""
    load(file::JSON{1})

Open a JSON file and parse its data in a dictionary.
"""
function FilesIO.load(file::JSON{1})
    open(FilesIO.filepath(file), "r") do f
        data = JSON3.read(f)
        return Dict(data)
    end
end

"""
    load(file::TXT{1})

Open a TEXT file and parse its data in a list of strings.
"""
function FilesIO.load(file::TXT{1})
    return readlines(FilesIO.filepath(file))
end

"""
    load(file::YAML{1})

Open a YAML file and parse its data in a dictionary.
"""
function FilesIO.load(file::YAML{1})
    return YAMLLib.load_file(FilesIO.filepath(file); dicttype=Dict{Symbol,Any})
end

end
