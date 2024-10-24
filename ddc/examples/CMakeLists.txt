set(ddc_EXAMPLES 
    2d_diffusive_convection
    3d_diffusive_convection
    2d_finger_convection
    3d_finger_convection
)

foreach (program ${ddc_EXAMPLES})
    install_channelflow_application(${program} examples)
    target_link_libraries(${program}_app PUBLIC ddc)
endforeach (program)
