target_sources(cfd
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Convection.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/HLLC.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Roe.cpp
)
