INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/lib/utgraphicsutil)
AUX_SOURCE_DIRECTORY(${CMAKE_SOURCE_DIR}/lib/utgraphicsutil libutgu_src)
FIND_PACKAGE(JPEG REQUIRED)
FIND_PACKAGE(PNG REQUIRED)
ADD_LIBRARY(utgraphicsutil STATIC ${libutgu_src})
TARGET_LINK_LIBRARIES(utgraphicsutil ${JPEG_LIBRARIES})
message("JPEG ${JPEG_INCLUDE_DIR}")
TARGET_INCLUDE_DIRECTORIES(utgraphicsutil SYSTEM BEFORE PRIVATE ${JPEG_INCLUDE_DIR})
TARGET_LINK_LIBRARIES(utgraphicsutil ${PNG_LIBRARIES})
message("JPEG ${PNG_INCLUDE_DIR}")
TARGET_INCLUDE_DIRECTORIES(utgraphicsutil SYSTEM BEFORE PRIVATE ${PNG_INCLUDE_DIR})
list(APPEND stdgl_libraries utgraphicsutil)