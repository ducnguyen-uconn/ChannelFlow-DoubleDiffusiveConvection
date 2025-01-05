#  
# Include and configure the channelflow tools
#
# This file is a part of channelflow version 2.0, https://channelflow.ch .
# License is GNU GPL version 2 or later: ./LICENSE
#
set(tool_APPS
  ddc_initialfield
  ddc_randomfield
  ddc_readEQ
  ddc_changegrid
  )

foreach (program ${tool_APPS})
    install_channelflow_application(${program} bin)
endforeach ()
