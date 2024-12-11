cmake_minimum_required(VERSION 3.27)
project(Stix C)
cmake_policy(SET CMP0135 NEW)
include(ExternalProject)


###############################################################################
# ## HSTLIB
# # build static library which is used by both giggle and stix
# set(HTSLIB_VERSION "1.21")
# set(HTS_SRC_DIR ${CMAKE_CURRENT_BINARY_DIR}/htslib)
# set(HTS_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/htslib)
# set(HTS_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/htslib)
#
# ExternalProject_Add(htslib
#     URL https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2
#     SOURCE_DIR ${HTS_SRC_DIR}
#     BINARY_DIR ${HTS_BINARY_DIR}
#     CONFIGURE_COMMAND COMMAND autoreconf -i 
#                       COMMAND ./configure --disable-bz2 --disable-lzma --disable-libcurl
#     BUILD_COMMAND make -C ${CMAKE_CURRENT_BINARY_DIR}/htslib
#     LOG_DOWNLOAD ON
#     LOG_CONFIGURE ON
#     INSTALL_COMMAND ""
#     PREFIX ""
# )
#
# add_library(hts STATIC IMPORTED)
# set_target_properties(hts PROPERTIES
#     IMPORTED_LOCATION ${HTS_BINARY_DIR}/libhts.a
#     INTERFACE_INCLUDE_DIRECTORIES ${HTS_INCLUDE_DIR}
# )
# add_dependencies(htslib hts)


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
ExternalProject_Add(sqlite_src
    URL http://www.sqlite.org/2017/sqlite-amalgamation-3170000.zip
    SOURCE_DIR ${CMAKE_BINARY_DIR}/sqlite
    PREFIX ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
)
add_library(sqlite STATIC ${CMAKE_CURRENT_BINARY_DIR}/sqlite/sqlite3.c)
set_target_properties(sqlite PROPERTIES
    IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/sqlite/libsqlite.a
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_BINARY_DIR}/sqlite
)

# ###############################################################################
## STIX
# The main event
set(HTS_INCLUDE_DIR ${CMAKE_BINARY_DIR}/giggle/htslib)
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













