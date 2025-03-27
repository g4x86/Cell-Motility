@ECHO OFF
REM Compile ProcSignalHandler using CMake and Visual Studio

SETLOCAL

REM Configure CMake C++ project by generating Visual Studio project files.
set CMAKE_EXEC=D:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe
"%CMAKE_EXEC%" -G "Visual Studio 17 2022" -A x64  -DCMAKE_CONFIGURATION_TYPES:STRING="Debug" -DCMAKE_INSTALL_PREFIX:PATH="C:\Users\granville\Documents\Repos\programs\CProjects\Cell-Motility\out\install\x64-Debug" -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON "C:\Users\granville\Documents\Repos\programs\CProjects\Cell-Motility"

REM Compile CMake C++ project
"%CMAKE_EXEC%" --build .

ENDLOCAL

@ECHO ON
