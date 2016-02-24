# RBSP

The following should be added to `~/.bashrc` and (?) `~/.bash_profile`

    function spedas {
      module load idl
      export ROOT_DATA_DIR="/export/scratch/users/mceachern/rbsp"
      ROOTDIR=~/Desktop/rbsp/packages/
      idl $* -IDL_PATH "+$ROOTDIR:<IDL_DEFAULT>"
    }
    export -f spedas

Bleeding-edge SPEDAS software comes from `http://themis.ssl.berkeley.edu/socware/bleeding_edge/`. 

Geopack library is located at `http://ampere.jhuapl.edu/code/idl_geopack.html`. Note that version 9.3 has dependencies that the physics department machines can't satisfy, but 7.6 seems to work. 

Icy, the IDL SPICE toolkit, comes from `http://naif.jpl.nasa.gov/naif/toolkit_IDL.html`. 

