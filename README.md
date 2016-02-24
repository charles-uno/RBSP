# RBSP

The following should be added to `~/.bashrc` and (?) `~/.bash_profile`

    function spedas {
      module load idl
      export ROOT_DATA_DIR="/export/scratch/users/mceachern/rbsp"
      ROOTDIR=~/Desktop/rbsp/packages/
      idl $* -IDL_PATH "+$ROOTDIR:<IDL_DEFAULT>"
    }
    export -f spedas
