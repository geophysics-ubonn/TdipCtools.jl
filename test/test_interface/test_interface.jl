@testset "convert_files" begin
    convert_files(
        "test_interface/t.txt",
        "test_interface/d.txt",
        "test_interface/cinv.txt",
        "test_interface/freqs.txt",
        "test_interface/spectrum.txt"
    )
end