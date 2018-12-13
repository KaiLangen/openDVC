add_library (colourizer
             colourizer.cpp
             mvMultiFrame.cpp
             mvSameColour.cpp
             mvSearchColour.cpp)

include (FILES colourizer.h mvMultiFrame.h mvSameColour.h mvSearchColour.h
         DESTINATION include)

install (TARGETS colourizer DESTINATION lib)

add_executable(colour main.cpp)

target_link_libraries(colour common colourizer)

                     
