set(
    ddc_APPS
    ddc_simulateflow
    ddc_findsoln
    ddc_continuesoln
    ddc_findeigenvals
    # ddc_edgetracking
    # ddc_addbaseflow
)

foreach (program ${ddc_APPS})
    install_channelflow_application(${program} bin)
    target_link_libraries(${program}_app PUBLIC ddc)
endforeach (program)

