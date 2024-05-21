

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



module load bioinfo-tools Nextflow/22.10.1

Note that NXF_HOME is set to $HOME/.nextflow
Please change NXF_HOME to a place in your project directory (export NXF_HOME=yourprojectfolder)

andbou@rackham3: (nfcore_mag) nf-core-mag-test: export NXF_HOME=/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/.nextflow
andbou@rackham3: (nfcore_mag) nf-core-mag-test: nextflow run nf-core/mag -profile test,uppmax --project naiss2024-5-1
N E X T F L O W  ~  version 22.10.2
Pulling nf-core/mag ...
 downloaded from https://github.com/nf-core/mag.git
Launching `https://github.com/nf-core/mag` [cheeky_pasteur] DSL2 - revision: a50117e4dd [master]
Nextflow version 22.10.2 does not match workflow required version: >=23.04.0


/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/databases/busco-db/bacteria_odb10


WARN: Singularity cache directory has not been defined -- Remote image will be stored in the path: /crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/nf-core-mag-test/work/singularity -- Use the environment variable NXF_SINGULARITY_CACHEDIR to specify a different location



nextflow run nf-core/mag -r 3.0.0 -profile uppmax --project naiss2024-5-1 --outdir CHST

nextflow run <pipeline> -bg > my-file.log