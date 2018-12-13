set (dec_src_files
     cavlcDec.cpp
     cavlcDec.h
     corrModel.cpp
     corrModel.h
     decoder.cpp
     decoder.h
     ldpcaDec.cpp
     ldpcaDec.h
     sideInformation.cpp
     sideInformation.h)

add_library(decoder ${dec_src_files})

include_directories(${CMAKE_SOURCE_DIR}/src/common)
target_link_libraries (decoder common)

add_executable(deDVC main.cpp)
target_link_libraries(deDVC common decoder)
