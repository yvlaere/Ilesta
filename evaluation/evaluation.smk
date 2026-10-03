import os
import random
import gzip

configfile: "config.yaml"

COVERAGES = config["coverages"]
prefix = config["output_prefix"]

# -------------------------------------------------------------------
# Rule all
# -------------------------------------------------------------------
rule all:
    input:
        "evaluation/output/n50_vs_coverage.png",
        "evaluation/output/largest_contig_vs_coverage.png",
        expand("evaluation/output/coverage_{cov}/subsample_stats.txt", cov=COVERAGES)


# -------------------------------------------------------------------
# Subsample reads
# -------------------------------------------------------------------
rule subsample_reads:
    input:
        fq=config["reads_fq"]
    output:
        fq="evaluation/output/subsampled_{cov}.fq"
    params:
        genome_size=config["genome_size"]
    run:

        def open_reads(path):
            if str(path).endswith(".gz"):
                return gzip.open(path, "rt")
            return open(path, "r")

        cov = float(wildcards.cov)
        rng = random.Random()

        target_bp = params.genome_size * cov

        # first pass: total bp
        total_bp = 0
        with open_reads(input.fq) as fin:
            while True:
                h = fin.readline()
                if not h:
                    break
                s = fin.readline()
                sep = fin.readline()
                q = fin.readline()
                if not q:
                    break
                total_bp += len(s.strip())

        keep_p = min(1.0, target_bp / total_bp)

        # second pass
        with open_reads(input.fq) as fin, open(output.fq, "w") as fout:
            while True:
                h = fin.readline()
                if not h:
                    break
                s = fin.readline()
                sep = fin.readline()
                q = fin.readline()
                if not q:
                    break

                if rng.random() < keep_p:
                    fout.write(h)
                    fout.write(s)
                    fout.write(sep)
                    fout.write(q)

rule read_subsample_stats:
    input:
        fq="evaluation/output/subsampled_{cov}.fq"
    output:
        stats="evaluation/output/coverage_{cov}/subsample_stats.txt"
    run:
        total_reads = 0
        total_bases = 0
        with open(input.fq) as f:
            while True:
                h = f.readline()
                if not h:
                    break
                s = f.readline()
                sep = f.readline()
                q = f.readline()
                if not q:
                    break
                total_reads += 1
                total_bases += len(s.strip())

        # read length distribution
        read_lengths = []
        with open(input.fq) as f:
            while True:
                h = f.readline()
                if not h:
                    break
                s = f.readline()
                sep = f.readline()
                q = f.readline()
                if not q:
                    break
                read_lengths.append(len(s.strip()))

        os.makedirs(os.path.dirname(output.stats), exist_ok=True)
        with open(output.stats, "w") as f:
            f.write(f"Total reads: {total_reads}\n")
            f.write(f"Total bases: {total_bases}\n")
            f.write(f"Read lengths: {read_lengths}\n")

# -------------------------------------------------------------------
# Assemble with Ilesta
# -------------------------------------------------------------------
rule assemble:
    input:
        fq="evaluation/output/subsampled_{cov}.fq"
    output:
        fa=f"evaluation/output/coverage_{{cov}}/{prefix}.fa"
    params:
        completion_flag=(
            "--completion-enabled=true"
            if config["completion_enabled"]
            else "--completion-enabled"
        )
    threads: config["threads"]
    shell:
        """
        mkdir -p evaluation/output/coverage_{wildcards.cov}

        {config[ilesta_path]} assemble \
            --reads-fq={input.fq} \
            --output-dir=evaluation/output/coverage_{wildcards.cov} \
            --threads={threads} \
            --paf={config[paf]} \
            --min-read-length={config[min_read_length]} \
            --min-base-quality={config[min_base_quality]} \
            --genome-size={config[genome_size]} \
            --min-overlap-length={config[min_overlap_length]} \
            --min-overlap-count={config[min_overlap_count]} \
            --min-percent-identity={config[min_percent_identity]} \
            --overhang-ratio={config[overhang_ratio]} \
            --output-prefix={config[output_prefix]} \
            --max-bubble-length={config[max_bubble_length]} \
            --min-support-ratio={config[min_support_ratio]} \
            --max-tip-len={config[max_tip_len]} \
            --fuzz={config[fuzz]} \
            --cleanup-iterations={config[cleanup_iterations]} \
            --short-edge-ratio={config[short_edge_ratio]} \
            {params.completion_flag} \
            --completion-rounds={config[completion_rounds]} \
            --completion-min-alignment-len={config[completion_min_alignment_len]} \
            --completion-min-identity={config[completion_min_identity]}
        """


# -------------------------------------------------------------------
# Compute assembly stats
# -------------------------------------------------------------------
rule compute_stats:
    input:
        fa=f"evaluation/output/coverage_{{cov}}/{prefix}.fa"
    output:
        stats="evaluation/output/coverage_{cov}/stats.txt"
    run:
        def get_lengths(fasta):
            lengths = []
            cur = 0
            with open(fasta) as f:
                for line in f:
                    if line.startswith(">"):
                        if cur:
                            lengths.append(cur)
                        cur = 0
                    else:
                        cur += len(line.strip())
                if cur:
                    lengths.append(cur)
            return lengths

        lengths = get_lengths(input.fa)

        if not lengths:
            n50 = 0
            largest = 0
        else:
            lengths_sorted = sorted(lengths, reverse=True)
            total = sum(lengths_sorted)

            cum = 0
            n50 = 0
            for l in lengths_sorted:
                cum += l
                if cum >= total / 2:
                    n50 = l
                    break

            largest = max(lengths_sorted)

        with open(output.stats, "w") as f:
            f.write(f"{wildcards.cov}\t{n50}\t{largest}\n")


# -------------------------------------------------------------------
# Aggregate stats
# -------------------------------------------------------------------
rule aggregate_stats:
    input:
        expand("evaluation/output/coverage_{cov}/stats.txt", cov=COVERAGES)
    output:
        "evaluation/output/metrics.tsv"
    run:
        rows = []
        for cov in COVERAGES:
            with open(f"evaluation/output/coverage_{cov}/stats.txt") as f:
                cov_, n50, largest = f.read().strip().split("\t")
                rows.append((float(cov_), int(n50), int(largest)))

        with open(output[0], "w") as f:
            f.write("coverage\tn50\tlargest_contig\n")
            for r in sorted(rows):
                f.write(f"{r[0]}\t{r[1]}\t{r[2]}\n")


# -------------------------------------------------------------------
# Plotting
# -------------------------------------------------------------------
rule plot:
    input:
        "evaluation/output/metrics.tsv"
    output:
        n50="evaluation/output/n50_vs_coverage.png",
        largest="evaluation/output/largest_contig_vs_coverage.png"
    run:
        import matplotlib.pyplot as plt

        cov, n50, largest = [], [], []

        with open(input[0]) as f:
            next(f)
            for line in f:
                c, n, l = line.strip().split("\t")
                cov.append(float(c))
                n50.append(int(n))
                largest.append(int(l))

        # N50 plot
        plt.figure()
        plt.plot(cov, n50, marker="o")
        plt.xlabel("Coverage")
        plt.ylabel("N50")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(output.n50, dpi=300)
        plt.close()

        # Largest contig plot
        plt.figure()
        plt.plot(cov, largest, marker="o")
        plt.xlabel("Coverage")
        plt.ylabel("Largest contig")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(output.largest, dpi=300)