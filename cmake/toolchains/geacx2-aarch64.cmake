set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR aarch64)

set(CMAKE_C_COMPILER "/usr/bin/aarch64-linux-gnu-gcc" CACHE FILEPATH "" FORCE)
set(CMAKE_CXX_COMPILER "/usr/bin/aarch64-linux-gnu-g++" CACHE FILEPATH "" FORCE)
set(ENV{BUILD_ARCH} "aarch64")

# FindPythonExtra otherwise queries the x86_64 host interpreter during cross
# compilation and gives ARM64 ROSIDL extensions an unloadable x86_64 suffix.
set(PYTHON_SOABI "cpython-310-aarch64-linux-gnu" CACHE INTERNAL
        "Target Python extension ABI" FORCE)

list(PREPEND CMAKE_PREFIX_PATH "/usr/lib/aarch64-linux-gnu/cmake")
list(PREPEND CMAKE_INCLUDE_PATH "/usr/include/aarch64-linux-gnu" "/usr/include")
list(PREPEND CMAKE_LIBRARY_PATH "/usr/lib/aarch64-linux-gnu")

set(Unwind_INCLUDE_DIR "/usr/include/aarch64-linux-gnu" CACHE PATH "" FORCE)
set(Unwind_LIBRARY "/usr/lib/aarch64-linux-gnu/libunwind.so" CACHE FILEPATH "" FORCE)
set(Unwind_PLATFORM_LIBRARY "/usr/lib/aarch64-linux-gnu/libunwind-aarch64.so" CACHE FILEPATH "" FORCE)
set(Qt5_DIR "/usr/lib/aarch64-linux-gnu/cmake/Qt5" CACHE PATH "" FORCE)

include("/opt/Cross-Compilation/toolchain.cmake")

# The vendor toolchain replaces CMAKE_TOOLCHAIN_FILE while it is being loaded.
# Pin it back to this overlay so try_compile and later reconfigure steps keep
# the target compiler and multiarch search paths.
set(CMAKE_TOOLCHAIN_FILE "${CMAKE_CURRENT_LIST_FILE}" CACHE FILEPATH "" FORCE)
