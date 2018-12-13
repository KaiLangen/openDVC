project (common)

set (com_src_files
     bitstream.cpp
     bitstream.h
     cavlc.cpp
     cavlc.h
     codec.cpp
     codec.h
     defs.h
     ldpca.cpp
     ldpca.h
     ME.cpp
     ME.h
     regExp.h
     transform.cpp
     transform.h
     types.h
)

add_library(common ${com_src_files})
                     
