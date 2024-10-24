set(ddc_VALIDATIONS
    yang2021jfm_case3_2d
)

foreach (program ${ddc_VALIDATIONS})
    install_channelflow_application(${program} validations)
    target_link_libraries(${program}_app PUBLIC ddc)
endforeach (program)
