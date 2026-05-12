class ConfigUtils {

    static int calculateCpus(memory) {
        def mem
        if (memory instanceof nextflow.util.MemoryUnit) {
            mem = memory.toMega()
        } else {
            throw new IllegalArgumentException("Memory not formatted correctly. Please specify memory as MemoryUnit.")
        }
        // Google Cloud Platform (GCP) restricts a maximum of 6 GB (6144 MB) provided per-CPU
        // Compute # of CPUs based on this (1 CPU minimum required)
        def cpus = Math.max(1, Math.ceil(mem / 6144).intValue())

        // # of CPUs required to be even by GCP for tasks with CPUs > 1
        if (cpus > 1) {
            cpus = cpus + (cpus % 2)
        }

        return cpus
    }
}
