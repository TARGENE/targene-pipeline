FROM julia:1.6.3

ENV JULIA_DEPOT_PATH=/opt

COPY . /EstimationPipeline.jl 

WORKDIR /EstimationPipeline.jl 

RUN julia --project -e"using Pkg; Pkg.instantiate(); Pkg.test()"