Given a set of single-sample VCF files, it counts how many times each variant is called in the set. 
It reports this in the following format: 

```
chrom  start  stop  alts  count_with_filter_pass  count_with_other_filters
chr1  11021  11022  ('A',)  1  0
```

### Build Docker

Build the Docker image using the following command after replacing `[TAG]` with your preferred tag. 

```
docker build --platform linux/amd64 -t [TAG] .
```


### Inputs

The following is an example input to the workflow.

```json
{
    "VariantOccurrenceFrequency.sample_ids": [
        "SP0001643",
        "SP0001677",
        "SP0001710",
        "SP0001952",
        "SP0002342"
    ],
    "VariantOccurrenceFrequency.vcf_files": [
        "SP0001643.vcf.gz",
        "SP0001677.vcf.gz",
        "SP0001710.vcf.gz",
        "SP0001952.vcf.gz",
        "SP0002342.vcf.gz"
    ],
    "VariantOccurrenceFrequency.runtime_override_encode": {"docker": "pzm:latest"},
    "VariantOccurrenceFrequency.runtime_override_merge": {"docker": "pzm:latest"},
    "VariantOccurrenceFrequency.runtime_override_decode": {"docker": "pzm:latest"}
}

```
