#!/bin/bash
## Job Name
#SBATCH --job-name=DSS
## Allocation Definition
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=10-00:00:00
## Memory per node
#SBATCH --mem=100G
#SBATCH --signal=USR2
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yaaminiv@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/yaaminiv/Hawes/analyses/DSS

################################################################
# Set user-defined variables inside this box.

# Set container path
container_path="/usr/lusers/yaaminiv"

# Set container name
container="rstudio-4.0.2.yrv-v1.0.sif"
################################################################

## MAKE NO CHANGES BELOW THIS LINE

module load singularity



export PASSWORD=$(openssl rand -base64 15)
# get unused socket per https://unix.stackexchange.com/a/132524
# tiny race condition between the python & singularity commands
readonly PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')
cat 1>&2 <<END
1. SSH tunnel from your workstation using the following command:

   ssh -N -L 8787:${HOSTNAME}:${PORT} ${USER}@mox.hyak.uw.edu

   and point your web browser to http://localhost:8787

2. log in to RStudio Server using the following credentials:

   user: ${USER}
   password: ${PASSWORD}

When done using RStudio Server, terminate the job by:

1. Exit the RStudio Session ("power" button in the top right corner of the RStudio window)
2. Issue the following command on the login node:

      scancel -f ${SLURM_JOB_ID}
END


# User-installed R packages go into their home directory
if [ ! -e ${HOME}/.Renviron ]
then
  printf '\nNOTE: creating ~/.Renviron file\n\n'
  echo 'R_LIBS_USER=~/R/%p-library/%v' >> ${HOME}/.Renviron
fi

# This example bind mounts the /gscratch/scrubbed/${USER} directory on the host into the Singularity container.
# By default the only host file systems mounted within the container are $HOME, /tmp, /proc, /sys, and /dev.
singularity exec \
--bind="$TMPDIR/var/lib:/var/lib/rstudio-server" \
--bind="$TMPDIR/var/run:/var/run/rstudio-server" \
--bind="$TMPDIR/tmp:/tmp" \
--bind=/gscratch/scrubbed/${USER} \
${container_path}/${container} \
rserver --www-port ${PORT} --auth-none=0 --auth-pam-helper-path=pam-helper --auth-stay-signed-in-days=30 --auth-timeout-minutes=0
