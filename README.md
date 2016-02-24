# RBSP

The following should be added to `~/.bashrc` 

  function spedas {
    module load idl
    export ROOT_DATA_DIR="/export/scratch/users/mceachern/RBSP"
    ROOTDIR=~/Desktop/RBSP/packages/
    idl $* -IDL_PATH "+$ROOTDIR:<IDL_DEFAULT>"
  }
  export -f spedas
