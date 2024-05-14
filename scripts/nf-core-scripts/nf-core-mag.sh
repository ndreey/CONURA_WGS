

# Do this before loading
export SINGULARITY_DISABLE_CACHE=true
export APPTAINER_DISABLE_CACHE=true
unset NXF_SINGULARITY_CACHEDIR

module load bioinfo-tools Nextflow/22.10.1 nf-core-pipelines/latest
module load bioinfo-tools Nextflow/22.10.1

# When one loads the Nextflow module it automatically resets back to $HOME/user so always call this
export NXF_HOME=/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/.nextflow


#nextflow run $NF_CORE_PIPELINES/mag/2.3.2/workflow -profile test,uppmax --project naiss2024-5-1

nf-core/2.6