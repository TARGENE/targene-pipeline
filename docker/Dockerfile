FROM julia:1.6.3

# Install GCTA: for GRM computation
RUN apt update -y && apt install -y zip

RUN curl -o gcta64dir.zip https://yanglab.westlake.edu.cn/software/gcta/bin/gcta_v1.94.0Beta_linux_kernel_3_x86_64.zip

RUN unzip gcta64dir.zip

RUN mv gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static /usr/local/bin/gcta64

RUN chmod a+x /usr/local/bin/gcta64

# Setup julia

ENV JULIA_DEPOT_PATH=/opt

COPY . /EstimationPipeline.jl 

WORKDIR /EstimationPipeline.jl 

RUN julia --project -e"using Pkg; Pkg.instantiate(); Pkg.test()"