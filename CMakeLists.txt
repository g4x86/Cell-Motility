# CMake configuration

# Set CMake version
cmake_minimum_required(VERSION 3.8)

# Project name
project(CellMotility
    LANGUAGES CXX
    VERSION 1.0.0
)

# Enable verbose build output
set(CMAKE_VERBOSE_MAKEFILE ON)

# Generate compile_commands.json with all compilation commands
if(${CMAKE_VERSION} VERSION_LESS "3.18.0")
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
else()
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

# Source files
set(PROJECT_SOURCES
    ${PROJECT_SOURCE_DIR}/include/typedefs.hpp
    ${PROJECT_SOURCE_DIR}/include/Actin.hpp
    ${PROJECT_SOURCE_DIR}/include/ADF.hpp
    ${PROJECT_SOURCE_DIR}/include/algorithms.hpp
    ${PROJECT_SOURCE_DIR}/include/ARP23.hpp
    ${PROJECT_SOURCE_DIR}/include/BranchTree.hpp
    ${PROJECT_SOURCE_DIR}/include/CAP.hpp
    ${PROJECT_SOURCE_DIR}/include/CellStatisticsCalculator.hpp
    ${PROJECT_SOURCE_DIR}/include/constants.hpp
    ${PROJECT_SOURCE_DIR}/include/Coordinate.hpp
    ${PROJECT_SOURCE_DIR}/include/DiscreteEvent.hpp
    ${PROJECT_SOURCE_DIR}/include/DiscreteEventSimulator.hpp
    ${PROJECT_SOURCE_DIR}/include/Event.hpp
    ${PROJECT_SOURCE_DIR}/include/EventSimulator.hpp
    ${PROJECT_SOURCE_DIR}/include/FilamentBranch.hpp
    ${PROJECT_SOURCE_DIR}/include/FilamentReaction.hpp
    ${PROJECT_SOURCE_DIR}/include/initializeParameterTable.hpp
    ${PROJECT_SOURCE_DIR}/include/initializeReactionTypeTable.hpp
    ${PROJECT_SOURCE_DIR}/include/InputFile.hpp
    ${PROJECT_SOURCE_DIR}/include/Line.hpp
    ${PROJECT_SOURCE_DIR}/include/MolecularDistributionFunction.hpp
    ${PROJECT_SOURCE_DIR}/include/MolecularDistribution.hpp
    ${PROJECT_SOURCE_DIR}/include/Molecule.hpp
    ${PROJECT_SOURCE_DIR}/include/MotileCell.hpp
    ${PROJECT_SOURCE_DIR}/include/OutputFile.hpp
    ${PROJECT_SOURCE_DIR}/include/ParameterTable.hpp
    ${PROJECT_SOURCE_DIR}/include/ReactionTypeTable.hpp
    ${PROJECT_SOURCE_DIR}/include/SpatialBoundary.hpp
    ${PROJECT_SOURCE_DIR}/include/SurfaceTopology.hpp
    ${PROJECT_SOURCE_DIR}/include/TokenIterator.hpp
    ${PROJECT_SOURCE_DIR}/include/Triangle.hpp
    ${PROJECT_SOURCE_DIR}/include/UniformMolecularDistributionFunction.hpp
    ${PROJECT_SOURCE_DIR}/include/UniformMolecularDistribution.hpp
    ${PROJECT_SOURCE_DIR}/include/Vector.hpp
    ${PROJECT_SOURCE_DIR}/include/VertexEdgeFacet.hpp
    ${PROJECT_SOURCE_DIR}/src/Actin.cpp
    ${PROJECT_SOURCE_DIR}/src/ADF.cpp
    ${PROJECT_SOURCE_DIR}/src/algorithms.cpp
    ${PROJECT_SOURCE_DIR}/src/ARP23.cpp
    ${PROJECT_SOURCE_DIR}/src/BranchTree.cpp
    ${PROJECT_SOURCE_DIR}/src/CAP.cpp
    ${PROJECT_SOURCE_DIR}/src/CellStatisticsCalculator.cpp
    ${PROJECT_SOURCE_DIR}/src/Coordinate.cpp
    ${PROJECT_SOURCE_DIR}/src/DiscreteEvent.cpp
    ${PROJECT_SOURCE_DIR}/src/DiscreteEventSimulator.cpp
    ${PROJECT_SOURCE_DIR}/src/Event.cpp
    ${PROJECT_SOURCE_DIR}/src/EventSimulator.cpp
    ${PROJECT_SOURCE_DIR}/src/FilamentBranch.cpp
    ${PROJECT_SOURCE_DIR}/src/FilamentReaction.cpp
    ${PROJECT_SOURCE_DIR}/src/initializeParameterTable.cpp
    ${PROJECT_SOURCE_DIR}/src/initializeReactionTypeTable.cpp
    ${PROJECT_SOURCE_DIR}/src/InputFile.cpp
    ${PROJECT_SOURCE_DIR}/src/Line.cpp
    ${PROJECT_SOURCE_DIR}/src/MolecularDistribution.cpp
    ${PROJECT_SOURCE_DIR}/src/MolecularDistributionFunction.cpp
    ${PROJECT_SOURCE_DIR}/src/Molecule.cpp
    ${PROJECT_SOURCE_DIR}/src/MotileCell.cpp
    ${PROJECT_SOURCE_DIR}/src/OutputFile.cpp
    ${PROJECT_SOURCE_DIR}/src/ParameterTable.cpp
    ${PROJECT_SOURCE_DIR}/src/ReactionTypeTable.cpp
    ${PROJECT_SOURCE_DIR}/src/SpatialBoundary.cpp
    ${PROJECT_SOURCE_DIR}/src/SurfaceTopology.cpp
    ${PROJECT_SOURCE_DIR}/src/TokenIterator.cpp
    ${PROJECT_SOURCE_DIR}/src/Triangle.cpp
    ${PROJECT_SOURCE_DIR}/src/UniformMolecularDistribution.cpp
    ${PROJECT_SOURCE_DIR}/src/UniformMolecularDistributionFunction.cpp
    ${PROJECT_SOURCE_DIR}/src/Vector.cpp
    ${PROJECT_SOURCE_DIR}/src/VertexEdgeFacet.cpp
    ${PROJECT_SOURCE_DIR}/src/main.cpp
)

# Main program
add_executable(${PROJECT_NAME}
    ${PROJECT_SOURCES}
)

# Set C++ standard version to c++98 in order to support
# a wide range of compilers and systems.
set_target_properties(${PROJECT_NAME}
    PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
)

# Add include directory
target_include_directories(${PROJECT_NAME}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)

# Set C++ language features for >= c++11
target_compile_features(${PROJECT_NAME}
    PUBLIC
        cxx_constexpr
        cxx_defaulted_functions
)

# Set compilation and link options
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # Set compilation options
    target_compile_options(${PROJECT_NAME}
        PRIVATE
            # Enable exception handling
            -fexceptions
            # Enable warnings for strict compliance with ISO C/C++
            -Wpedantic
            # Enable most warnings
            -Wall
            # Enable extra warnings
            -Wextra
            # Turn all warnings into errors
            # -Werror
    )
    # Set link options for GCC/Clang in Unix-like systems
    if(NOT Win32)
        if(APPLE)
            if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
                target_link_options(${PROJECT_NAME}
                    PRIVATE
                    -stdlib=libc++
                )
            elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
                target_link_options(${PROJECT_NAME}
                    PRIVATE
                    # libgcc is only available for GCC in macOS
                    -static-libgcc
                )
            endif()
        elseif(UNIX)
            target_link_options(${PROJECT_NAME}
                PRIVATE
                # libgcc is available in Linux?Unix regardless of compilers
                -static-libgcc
            )
        endif()
        target_link_options(${PROJECT_NAME}
            PRIVATE
            -static-libstdc++
        )
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # Set compilation options
    target_compile_options(${PROJECT_NAME}
        PRIVATE
            # Standard C++ exception handling
            # Ref: https://docs.microsoft.com/en-us/cpp/build/reference/eh-exception-handling-model?view=msvc-170
            /EHsc
            # Enable warnings for strict compliance with ISO C/C++
            /permissive-
            # Enable __cplusplus macro
            /Zc:__cplusplus
            # Enable all warnings including /W4 and those off by default
            /Wall
            # Turn all warnings into errors
            # /WX
    )
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "OpenWatcom")
    # Set compilation options
    target_compile_options(${PROJECT_NAME}
        PRIVATE
            # Enable exception handling
            -xs
            # Enable warnings for strict compliance with ISO C/C++
            -za
            # Enable all warnings
            -wx
            # Turn all warnings into errors
            # -we
    )
endif()
