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
            runtime_override = runtime_override_decode
    }

    output {
        File variants_frequence = decode.variants_frequence
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
        disks: 100,
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
                    key = base64.b64encode(f"{variant.chrom}:{str(variant.start)}:{str(variant.stop)}:{str(variant.alts)}".encode("utf-8")).decode("utf-8")
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
        disks: 100,
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
                out_file.write(f"{k}\t{v[0]}\t{v[1]}\n")
        CODE
    >>>
}

task DecodeVariants {
    input {
        File encoded_variants_dict
        RuntimeAttr runtime_override
    }

    output {
        File variants_frequence = "variants.csv.gz"
    }

    RuntimeAttr runtime_default = object {
        cpu: 1,
        memory: 64,
        disks: 100,
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

        with gzip.open("variants.csv.gz", "wt", compresslevel=4) as out_file:
            out_file.write("\t".join(["chrom", "start", "stop", "alts", "count-non-pass-filter", 'count-pass-filter']))
            for variant, frequency in variants.items():
                x = base64.b64decode(variant).decode("utf-8")
                x = x.split(":")
                x.extend(frequency)
                out_file.write("\t".join([str(c) for c in x]) + "\n")
        CODE
    >>>
}
