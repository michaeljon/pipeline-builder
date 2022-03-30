# Ovation virome processing pipeline

## Human genome processing

(Note, this repo also contains our human genome pipeline, but it's been on the back burner for several weeks as we've been chasing down sars-cov2 variant data.)


## SARS-Cov-2 virome processing

The shell script `sars-cov2-pe-runner` is a wrapper around the python-based pipeline builder. It can be used for building and executing a single run or for passing into GNU parallel (`parallel ./sars-cov2-pe-runner < ~/samples_to_run`). 

Both python scripts will output their command line options. Most of them are optional with the exception of `sample` and `work-dir`. Input files are assumed to use a `<sample>_R[12].fastq.gz` format. The pipeline assumes use of our of our pipeline AMIs. These have the current tools in `~/bin` and the current reference data in `~/reference`.