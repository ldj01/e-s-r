install(PROGRAMS updatelads.py DESTINATION bin COMPONENT lasrc-aux
        EXCLUDE_FROM_ALL)
install(FILES api_interface.py config_utils.py DESTINATION bin
        COMPONENT lasrc-aux EXCLUDE_FROM_ALL)
