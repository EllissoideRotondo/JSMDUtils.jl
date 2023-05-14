
@testset "Format" verbose = true begin

    # Test CamelCase formatting
    @test format_camelcase("hello world") == "HelloWorld"
    @test format_camelcase("hello_world") == "HelloWorld"
    @test format_camelcase("hello") == "Hello"
    @test format_camelcase("test this_string") == "TestThisString"
    @test format_camelcase("Test   This") == "TestThis"

    # Test SnakeCase formatting 
    @test format_snakecase("hello world") == "hello_world"
    @test format_snakecase("HelloWorld") == "hello_world"
    @test format_snakecase("Hello_world") == "hello_world"
end
