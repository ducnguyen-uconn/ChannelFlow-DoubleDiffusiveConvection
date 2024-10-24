set(ddc_TESTS ddc_timeIntegrationTest)

file(
    COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/data/uinit.nc ${CMAKE_CURRENT_SOURCE_DIR}/data/ufinal.nc
    DESTINATION
    ${CMAKE_CURRENT_BINARY_DIR}/data
)

file(
    COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/data/tinit.nc ${CMAKE_CURRENT_SOURCE_DIR}/data/tfinal.nc
    DESTINATION
    ${CMAKE_CURRENT_BINARY_DIR}/data
)

foreach (program ${ddc_TESTS})
    install_channelflow_application(${program} OFF)
    target_link_libraries(${program}_app PUBLIC ddc)
endforeach (program)

add_serial_test(ddc_timeIntegration ddc_timeIntegrationTest)

if (USE_MPI)
    add_mpi_test(mpi_ddc_timeIntegration ddc_timeIntegrationTest)
endif ()
