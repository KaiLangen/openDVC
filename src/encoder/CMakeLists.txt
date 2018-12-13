set (enc_src_files
     encoder3channel.cpp
     encoder.h
     cavlcEnc.cpp
     cavlcEnc.h
     ldpcaEnc.cpp
     ldpcaEnc.h)

add_library (encoder ${enc_src_files})

include_directories (${CMAKE_SOURCE_DIR}/src/common)
target_link_libraries (encoder common)

add_executable(enDVC main.cpp)
target_link_libraries(enDVC common encoder)
