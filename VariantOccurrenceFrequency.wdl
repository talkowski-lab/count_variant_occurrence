version 1.0

import "GetShardInputs.wdl"

struct RuntimeAttr {
    Int? cpu
    Float? memory
    Int? disks
    Int? bootDiskSizeGb
    Int? preemptible
    Int? maxRetries
    String docker
}

workflow VariantOccurrenceFrequency {
    input {
        Array[String] sample_ids
        Array[File] vcf_files

        Int samples_per_shard = 2

        RuntimeAttr runtime_override_encode
        RuntimeAttr runtime_override_merge
        RuntimeAttr runtime_override_decode
        RuntimeAttr runtime_override_sortindex
    }

    Array[Pair[String, File]] zipped = zip(sample_ids, vcf_files)
    Int vcfs_length = length(zipped)

    Int num_samples = length(vcf_files)
    Float num_samples_float = num_samples
    Int num_shards = ceil(num_samples_float / samples_per_shard)

    scatter (i in range(num_shards)) {
        call GetShardInputs.GetShardInputs as GetShardVCFs {
            input:
                items_per_shard = samples_per_shard,
                shard_number = i,
                num_items = num_samples,
                all_items = vcf_files
        }

        call GetShardInputs.GetShardInputs as GetShardIds {
            input:
                items_per_shard = samples_per_shard,
                shard_number = i,
                num_items = num_samples,
                all_items = sample_ids
        }

        call EncodeVariants as encode_variants {
            input:
                vcf_files = GetShardVCFs.shard_items,
                sample_ids = GetShardIds.shard_items,
                runtime_override = runtime_override_encode
        }
    }

    call Merge as merge {
        input:
            encoded_files = encode_variants.encoded_files,
            runtime_override = runtime_override_merge
    }

    call DecodeVariants as decode {
        input:
            encoded_variants_dict = merge.encoded_variants_dict,
            input_vcfs_count = num_samples,
            runtime_override = runtime_override_decode
    }

    call SortCompressIndex as sortCompressIndex {
        input:
            variants_frequence = decode.variants_frequence,
            runtime_override = runtime_override_sortindex
    }

    output {
        File variants_frequency = sortCompressIndex.variants_frequency
        File variants_frequency_index = sortCompressIndex.variants_frequency_index
    }
}

task EncodeVariants {
    input {
        Array[File] vcf_files
        Array[String] sample_ids
        RuntimeAttr runtime_override
    }

    output {
        Array[File] encoded_files = glob("*.csv.gz")
    }

    RuntimeAttr runtime_default = object {
        cpu: 1,
        memory: 3.75,
        disks: 25 + (ceil(size(vcf_files, "GiB")) * 2),
        bootDiskSizeGb: 10,
        preemptible: 3,
        maxRetries: 1,
        docker: "python:slim-buster"
    }
    RuntimeAttr runtime_attr = select_first([runtime_override, runtime_default])

    runtime {
        cpu: select_first([runtime_attr.cpu, runtime_default.cpu])
        memory: select_first([runtime_attr.memory, runtime_default.memory]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disks, runtime_default.disks]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.bootDiskSizeGb, runtime_default.bootDiskSizeGb])
        preemptible: select_first([runtime_attr.preemptible, runtime_default.preemptible])
        maxRetries: select_first([runtime_attr.maxRetries, runtime_default.maxRetries])
        docker: select_first([runtime_attr.docker, runtime_default.docker])
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
        import base64
        import gzip
        import json
        import pysam

        vcf_files = ["~{sep='", "' vcf_files}"]
        sample_ids = ["~{sep='", "' sample_ids}"]
        for sample_id, filename in zip(sample_ids, vcf_files):
            vcf = pysam.VariantFile(filename)
            with gzip.open(f"{sample_id}.csv.gz", "wt", compresslevel=4) as out_file:
                for variant in vcf:
                    if "multiallelic" in variant.filter.keys():
                        pass
                    key = base64.b64encode(f"{variant.chrom.removeprefix('chr')}:{str(variant.pos)}:{variant.ref}:{variant.alts[0]}".encode("utf-8")).decode("utf-8")
                    filter_key = "1" if "PASS" in variant.filter.keys() else "0"
                    out_file.write(f"{key}\t{filter_key}\n")
        CODE
    >>>
}

task Merge {
    input {
        Array[Array[File]] encoded_files
        RuntimeAttr runtime_override
    }

    output {
        File encoded_variants_dict = "dict.csv.gz"
    }

    RuntimeAttr runtime_default = object {
        cpu: 1,
        memory: 64,
        disks: 25 + (ceil(size(flatten(encoded_files), "GiB")) * 2),
        bootDiskSizeGb: 10,
        preemptible: 3,
        maxRetries: 1,
        docker: "python:slim-buster"
    }
    RuntimeAttr runtime_attr = select_first([runtime_override, runtime_default])

    runtime {
        cpu: select_first([runtime_attr.cpu, runtime_default.cpu])
        memory: select_first([runtime_attr.memory, runtime_default.memory]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disks, runtime_default.disks]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.bootDiskSizeGb, runtime_default.bootDiskSizeGb])
        preemptible: select_first([runtime_attr.preemptible, runtime_default.preemptible])
        maxRetries: select_first([runtime_attr.maxRetries, runtime_default.maxRetries])
        docker: select_first([runtime_attr.docker, runtime_default.docker])
    }

    command <<<
        set -euo pipefail

        CSVS="~{write_json(encoded_files)}"

        python3 <<CODE
        import gzip
        import json
        from collections import defaultdict

        with open("$CSVS", "r") as file:
            data = json.load(file)
        filenames = [x for xs in data for x in xs]

        variants = defaultdict(lambda: [0, 0])
        for filename in filenames:
            with gzip.open(filename, "rt") as file:
                for line in file:
                    cols = line.strip().split()
                    variants[cols[0]][int(cols[1])] += 1

        with gzip.open("dict.csv.gz", "wt", compresslevel=4) as out_file:
            for k, v in variants.items():
                out_file.write(f"{k}\t{v[0] + v[1]}\t{v[1]}\n")
        CODE
    >>>
}

task DecodeVariants {
    input {
        File encoded_variants_dict
        Int input_vcfs_count
        RuntimeAttr runtime_override
    }

    output {
        File variants_frequence = "variants.tab.gz"
    }

    RuntimeAttr runtime_default = object {
        cpu: 1,
        memory: 64,
        disks: 25 + (ceil(size(encoded_variants_dict, "GiB")) * 2),
        bootDiskSizeGb: 10,
        preemptible: 3,
        maxRetries: 1,
        docker: "python:slim-buster"
    }
    RuntimeAttr runtime_attr = select_first([runtime_override, runtime_default])

    runtime {
        cpu: select_first([runtime_attr.cpu, runtime_default.cpu])
        memory: select_first([runtime_attr.memory, runtime_default.memory]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disks, runtime_default.disks]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.bootDiskSizeGb, runtime_default.bootDiskSizeGb])
        preemptible: select_first([runtime_attr.preemptible, runtime_default.preemptible])
        maxRetries: select_first([runtime_attr.maxRetries, runtime_default.maxRetries])
        docker: select_first([runtime_attr.docker, runtime_default.docker])
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
        import base64
        import gzip

        variants = {}
        with gzip.open("~{encoded_variants_dict}", "rt") as file:
            for line in file:
                line = line.strip().split("\t")
                variants[line[0]] = [line[1], line[2]]

        sample_count = int(~{input_vcfs_count})
        with gzip.open("variants.tab.gz", "wt", compresslevel=4) as out_file:
            out_file.write("\t".join(["#chrom", "pos", "ref", "alt", "all_cohort_count", "count_pass_filter", "sample_count", "all_cohort_af", "pass_cohort_af"]) + "\n")
            for variant, frequency in variants.items():
                x = base64.b64decode(variant).decode("utf-8")
                x = x.split(":")
                x.append(frequency[0])
                x.append(frequency[1])
                x.append(str(sample_count))
                x.append(float(int(frequency[0]) / sample_count * 100.0))  # all_cohort_af
                x.append(float(int(frequency[1]) / sample_count * 100.0))  # pass_cohort_af
                out_file.write(f"chr{x[0]}\t" + "\t".join([str(c) for c in x[1:]]) + "\n")
        CODE
    >>>
}

task SortCompressIndex {
    input {
        File variants_frequence
        RuntimeAttr runtime_override
    }

    output {
        File variants_frequency = "sorted_variants.tab.gz"
        File variants_frequency_index = "sorted_variants.tab.gz.tbi"
    }

    RuntimeAttr runtime_default = object {
        cpu: 1,
        memory: 64,
        disks: 100 + (ceil(size(variants_frequence, "GiB")) * 6),
        bootDiskSizeGb: 50,
        preemptible: 3,
        maxRetries: 1,
        docker: "python:slim-buster"
    }
    RuntimeAttr runtime_attr = select_first([runtime_override, runtime_default])

    runtime {
        cpu: select_first([runtime_attr.cpu, runtime_default.cpu])
        memory: select_first([runtime_attr.memory, runtime_default.memory]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disks, runtime_default.disks]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.bootDiskSizeGb, runtime_default.bootDiskSizeGb])
        preemptible: select_first([runtime_attr.preemptible, runtime_default.preemptible])
        maxRetries: select_first([runtime_attr.maxRetries, runtime_default.maxRetries])
        docker: select_first([runtime_attr.docker, runtime_default.docker])
    }

    command <<<
        set -euo pipefail

        gunzip -c ~{variants_frequence} > variants.tab
        (head -n 1 variants.tab && tail -n +2 variants.tab | sort -k1,1V -k2,2n) > sorted_variants.tab

        bgzip -c "sorted_variants.tab" > "sorted_variants.tab.gz"
        tabix -s 1 -b 2 -e 2 "sorted_variants.tab.gz"
    >>>
}
