module TestFromParamFile

using Test
using DataFrames
using TMLE
using Arrow
using JLD2
using TMLECLI
using YAML

# "local" profile assumes singularity is installed
args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"]

include("utils.jl")

@testset "Test ukb_estimands_file.config" begin
    cmd = `nextflow run main.nf -c test/configs/ukb_estimands_file.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0
    
    # Check QQ plot is here
    @test isfile(joinpath("results", "QQ.png"))

    # Check full HDF5 results
    results = jldopen(io -> io["results"], joinpath("results", "results.hdf5"))
    @test size(results, 1) > 40
    failed_results = retrieve_failed_results(results)
    ## All fails are due to fluctuation failure due non positive definite matrix
    ## This does not affect the OSE
    @test isempty(failed_results.OSE_GLM_GLM)
    ## Less than 1/3 affected: this is still quite significant
    @test size(failed_results.TMLE_GLM_GLM, 1) / size(results, 1) < 1/3
    @test all(startswith(x.msg, "Could not fluctuate") for x ∈ failed_results.TMLE_GLM_GLM)
    dataset = DataFrame(Arrow.Table(joinpath("results", "datasets", "all_genotypes.data.arrow")))
    check_fails_are_extremely_rare_traits(failed_results.TMLE_GLM_GLM, dataset; ncases=3)

    # Check summary file
    summary_results = YAML.load_file(joinpath("results", "results.summary.yaml"))
    @test size(summary_results, 1) == size(results, 1)

    # Check properly resumed
    resume_time = @elapsed run(cmd)
    @test resume_time < 1000
end

end

true