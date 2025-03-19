@testset "convert_files" begin
    convert_files(
        "test_interface/t.txt",
        "test_interface/d.txt",
        "test_interface/cinv.txt",
        "test_interface/freqs.txt",
        "test_interface/spectrum.txt",
        "test_interface/r0.txt"
    )

    convert_files(
        "test_interface/tomo_t.txt",
        "test_interface/tomo_d.txt",
        "test_interface/tomo_cinv.txt",
        "test_interface/tomo_freqs.txt",
        "test_interface/tomo_spectrum.txt",
        "test_interface/tomo_r0.txt"
    )
end