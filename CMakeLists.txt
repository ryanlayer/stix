cmake_minimum_required(VERSION 3.27)
project(Stix C)
cmake_policy(SET CMP0135 NEW)
include(ExternalProject)
include(FetchContent)


###############################################################################
## GIGGLE/HTSLIB
set(GIGGLE_BRANCH "cmake-build")
ExternalProject_Add(giggle_src
    GIT_REPOSITORY https://github.com/ryanlayer/giggle.git
    GIT_TAG ${GIGGLE_BRANCH}
    GIT_SUBMODULES ""
    PREFIX ""
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/giggle
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
)

add_library(libgiggle STATIC IMPORTED)
add_dependencies(libgiggle giggle_src)
set_target_properties(libgiggle PROPERTIES
    IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/giggle/lib/libgiggle.a
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_BINARY_DIR}/giggle/src
    OUTPUT_NAME giggle
)

add_library(hts STATIC IMPORTED)
set_target_properties(hts PROPERTIES
    IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/giggle/htslib/libhts.a
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_BINARY_DIR}/giggle/htslib
)

###############################################################################
## SQLITE
# ExternalProject_Add(sqlite_src
FetchContent_Declare(sqlite_src
    URL http://www.sqlite.org/2017/sqlite-amalgamation-3170000.zip
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/sqlite
    PREFIX ""
    # CONFIGURE_COMMAND ""
    # BUILD_COMMAND ""
    # INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
)
FetchContent_MakeAvailable(sqlite_src)
add_library(sqlite STATIC ${CMAKE_CURRENT_BINARY_DIR}/sqlite/sqlite3.c)
set_target_properties(sqlite PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_BINARY_DIR}/sqlite
)

# ###############################################################################
## STIX
# The main event
set(HTS_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/giggle/htslib)
set(GIGGLE_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/giggle/src)
# set(SQLITE_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/sqlite/)
set(SQLITE_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/)
set(STIX_SRC ${CMAKE_CURRENT_BINARY_DIR}/src)

set(STIX_SRC_FILES
    "${STIX_SRC}/ped.c"
    "${STIX_SRC}/search.c"
    "${STIX_SRC}/sharding_utils.c"
)
add_library(stix_obj OBJECT ${STIX_SRC_FILES})
target_include_directories(stix_obj PRIVATE
     ${HTS_INCLUDE_DIR}
     ${GIGGLE_INCLUDE_DIR}
     ${SQLITE_INCLUDE_DIR}
)

add_executable(stix ${STIX_SRC}/stix.c)
include_directories(stix PRIVATE
    ${HTS_INCLUDE_DIR}
    ${GIGGLE_INCLUDE_DIR}
    ${SQLITE_INCLUDE_DIR}
)
target_link_libraries(stix PRIVATE libgiggle hts sqlite dl z m pthread stix_obj)
set_target_properties(stix PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)













