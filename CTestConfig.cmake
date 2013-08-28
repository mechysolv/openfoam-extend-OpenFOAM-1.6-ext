## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
## # The following are required to uses Dart and the Cdash dashboard
##   ENABLE_TESTING()
##   INCLUDE(CTest)
set(CTEST_PROJECT_NAME "OpenFOAM-1.6-ext")
set(CTEST_NIGHTLY_START_TIME "00:00:00 EST")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "openfoam-extend.sourceforge.net")
set(CTEST_DROP_LOCATION "/CDash/submit.php?project=OpenFOAM-1.6-ext")
set(CTEST_DROP_SITE_CDASH TRUE)

## Run ctest in parallel if environment variable WM_NCOMPPROCS is set
IF (NOT $ENV{WM_NCOMPPROCS} STREQUAL "")
    # Will run ctest in parallel over $WM_NCOMPPROCS processors
    set(CMAKE_CTEST_COMMAND ${CMAKE_CTEST_COMMAND} --parallel $ENV{WM_NCOMPPROCS})
    MESSAGE("Running tests in parallel using $ENV{WM_NCOMPPROCS} processors")
ENDIF  (NOT $ENV{WM_NCOMPPROCS} STREQUAL "")

