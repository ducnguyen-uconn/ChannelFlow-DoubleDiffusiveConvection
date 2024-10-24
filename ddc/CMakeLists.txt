#
# 
# Include and configure the Double-Diffusive Convection module of channelflow.
#
# Original author: Duc Nguyen
#
set(
    ddc_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/ddcflags.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/dde.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ddc.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ddcdsi.cpp
)

set(
    ddc_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/ddcflags.h
    ${CMAKE_CURRENT_SOURCE_DIR}/dde.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ddc.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ddcdsi.h
)

# Define the target with appropriate dependencies
install_channelflow_library(ddc)
target_link_fftw(ddc PRIVATE)
target_link_libraries(ddc PUBLIC chflow)
if (WITH_NSOLVER)
    target_link_libraries(ddc PUBLIC nsolver)
endif ()

# Install header files
install(FILES ${ddc_HEADERS} DESTINATION include/modules/ddc)
