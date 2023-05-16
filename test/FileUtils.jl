
import JSMDUtils.FileUtils 
import JSMDInterfaces.FilesIO as jIO

@testset "IO" verbose=true begin 

    # Test TXT file
    path = "test/assets/file_txt.txt"
    fileTXT = FileUtils.TXT(path)

    @test jIO.filepath(fileTXT) == path
    
    data = jIO.load(fileTXT)
    @test data[1] == "Hello world. "
    @test data[2] == ""
    @test data[3] == "TXT file reading test."

    # Test JSON file
    path = "test/assets/file_test.json"
    fileJSON = FileUtils.JSON(path)

    @test jIO.filepath(fileJSON) == path 
    
    data = jIO.load(fileJSON)

    @test data[:filename][:title] == "json text"
    @test data[:filename][:number] == 1

    # Test YAML file 
    path = ["test/assets/file_test.yaml"]
    fileYAML = FileUtils.YAML(path...)

    @test jIO.filepath(fileYAML) == path[1]

    data = jIO.load(fileYAML)

    @test data[:name] == "TestFile1"
    @test data[:number] == 1
    @test data[:text] == ["hello", "world"]

end;
