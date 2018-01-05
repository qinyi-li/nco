:: cmake build all NCO dependencies obtained with 'git clone'
:: Pedro Vicente

@echo off
if not defined DevEnvDir (
 echo "%VS140COMNTOOLS%VsDevCmd.bat" 
 call "%VS140COMNTOOLS%VsDevCmd.bat" 
 echo "%VCINSTALLDIR%vcvarsall.bat" amd64
 call "%VCINSTALLDIR%vcvarsall.bat" amd64
 if errorlevel 1 goto :eof
)

if "%~1" == "crt" (
  set STATIC_CRT=ON
) else (
  set STATIC_CRT=OFF
)
set MSVC_VERSION="Visual Studio 14 2015 Win64"
echo using static crt %STATIC_CRT%
echo using %MSVC_VERSION%

:: replace the character string '\' with '/' needed for cmake
set root_win=%cd%
set root=%root_win:\=/%
echo cmake root is %root%

if exist %root_win%\zlib\build\zlib.sln (
 echo skipping zlib build
 goto build_hdf5
) else (
  echo building zlib
  pushd zlib
  mkdir build
  pushd build
  cmake .. -G %MSVC_VERSION% ^
           -DMSVC_USE_STATIC_CRT=%STATIC_CRT% ^
           -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON ^
           -DCMAKE_BUILD_TYPE=Debug ^
           -DBUILD_SHARED_LIBS=OFF
  msbuild zlib.sln /target:build /property:configuration=debug
  cp %root%\zlib\build\zconf.h %root%\zlib
  popd
  popd
  if errorlevel 1 goto :eof
)


:build_hdf5
if exist %root_win%\hdf5\build\bin\Debug\h5dump.exe (
 echo skipping hdf5 build
 goto build_curl
) else (
  echo building hdf5
  pushd hdf5
  mkdir build
  pushd build
  cmake .. -G %MSVC_VERSION% ^
           -DBUILD_STATIC_CRT_LIBS=%STATIC_CRT% ^
           -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON ^
           -DCMAKE_BUILD_TYPE=Debug ^
           -DBUILD_SHARED_LIBS=OFF ^
           -DBUILD_STATIC_EXECS=ON ^
           -DBUILD_TESTING=OFF ^
           -DHDF5_BUILD_EXAMPLES=OFF ^
           -DHDF5_BUILD_CPP_LIB=OFF ^
           -DHDF5_ENABLE_Z_LIB_SUPPORT=ON ^
           -DH5_ZLIB_HEADER=%root%/zlib/zlib.h ^
           -DZLIB_STATIC_LIBRARY:FILEPATH=%root%/zlib/build/Debug/zlibstaticd.lib ^
           -DZLIB_INCLUDE_DIRS:PATH=%root%/zlib
  msbuild HDF5.sln /target:build /property:configuration=debug
  popd
  popd
  if errorlevel 1 goto :eof
)


:build_curl
if exist %root_win%\curl\builds\libcurl-vc14-x64-debug-static-ipv6-sspi-winssl\lib\libcurl_a_debug.lib (
 echo skipping curl build
 goto build_netcdf
) else (
  echo building curl
  pushd curl
  call buildconf.bat
  pushd winbuild
  @echo on
  if %STATIC_CRT% == ON (
   nmake /f Makefile.vc mode=static vc=14 debug=yes gen_pdb=yes MACHINE=x64 RTLIBCFG=static
  ) else (
   nmake /f Makefile.vc mode=static vc=14 debug=yes gen_pdb=yes MACHINE=x64
  )
  @echo off
  popd
  popd
  if errorlevel 1 goto :eof
)


:build_netcdf
if exist %root_win%\netcdf-c\build\ncdump\ncdump.exe (
 echo skipping netcdf build
 goto test_netcdf
) else (
  echo building netcdf
  pushd netcdf-c
  mkdir build
  pushd build
  cmake .. -G %MSVC_VERSION% ^
           -DNC_USE_STATIC_CRT=%STATIC_CRT% ^
           -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON ^
           -DENABLE_TESTS=OFF ^
           -DCMAKE_BUILD_TYPE=Debug ^
           -DBUILD_SHARED_LIBS=OFF ^
           -DHDF5_HL_LIBRARY=%root%/hdf5/build/bin/Debug/libhdf5_hl_D.lib ^
           -DHDF5_C_LIBRARY=%root%/hdf5/build/bin/Debug/libhdf5_D.lib ^
           -DHDF5_INCLUDE_DIR=%root%/hdf5/src ^
           -DZLIB_LIBRARY:FILE=%root%/zlib/build/Debug/zlibstaticd.lib ^
           -DZLIB_INCLUDE_DIR:PATH=%root%/zlib ^
           -DHAVE_HDF5_H=%root%/hdf5/build ^
           -DHDF5_HL_INCLUDE_DIR=%root%/hdf5/hl/src ^
           -DCURL_LIBRARY=%root%/curl/builds/libcurl-vc14-x64-debug-static-ipv6-sspi-winssl/lib/libcurl_a_debug.lib ^
           -DCURL_INCLUDE_DIR=%root%/curl/include
  msbuild netcdf.sln /target:build /property:configuration=debug
  popd
  popd
  if errorlevel 1 goto :eof
)


:test_netcdf
if exist %root_win%\netcdf-c\build\ncdump\ncdump.exe (
 echo testing netcdf build
 %root_win%\netcdf-c\build\ncdump\ncdump.exe -h http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/cmap/enh/precip.mon.mean.nc
 goto build_expat
)

:build_expat
if exist %root_win%\libexpat\expat\build\expat.sln (
 echo skipping expat build
 goto build_udunits
) else (
  echo building expat
  pushd libexpat
  pushd expat
  mkdir build
  pushd build
  cmake .. -G %MSVC_VERSION% ^
           -DMSVC_USE_STATIC_CRT=%STATIC_CRT% ^
           -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON ^
           -DBUILD_SHARED_LIBS=OFF ^
           -DBUILD_shared=OFF
  msbuild expat.sln /target:build /property:configuration=debug
  popd
  popd
  popd
  if errorlevel 1 goto :eof
)


:build_udunits
if exist %root_win%\UDUNITS-2\build\udunits.sln (
 echo skipping udunits build
 goto build_nco
) else (
  echo building udunits
  pushd UDUNITS-2
  mkdir build
  pushd build
  cmake .. -G %MSVC_VERSION% ^
           -DMSVC_USE_STATIC_CRT=%STATIC_CRT% ^
           -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON ^
           -DBUILD_SHARED_LIBS=OFF ^
           -DEXPAT_INCLUDE_DIR=%root%/libexpat/expat/lib ^
           -DEXPAT_LIBRARY=%root%/libexpat/expat/build/Debug/expatd.lib
  msbuild udunits.sln /target:build /property:configuration=debug
  popd
  popd
  if errorlevel 1 goto :eof
)

:build_nco
if exist Debug\ncks.exe (
 echo skipping nco build
 goto test_nco
) else (
  rm -rf CMakeCache.txt CMakeFiles
  cmake .. -G %MSVC_VERSION% ^
  -DNCO_MSVC_USE_MT=%STATIC_CRT% ^
  -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON ^
  -DNETCDF_INCLUDE:PATH=%root%/netcdf-c/include ^
  -DNETCDF_LIBRARY:FILE=%root%/netcdf-c/build/liblib/Debug/netcdf.lib ^
  -DHDF5_LIBRARY:FILE=%root%/hdf5/build/bin/Debug/libhdf5_D.lib ^
  -DHDF5_HL_LIBRARY:FILE=%root%/hdf5/build/bin/Debug/libhdf5_hl_D.lib ^
  -DZLIB_LIBRARY:FILE=%root%/zlib/build/Debug/zlibstaticd.lib ^
  -DCURL_LIBRARY:FILE=%root%/curl/builds/libcurl-vc14-x64-debug-static-ipv6-sspi-winssl/lib/libcurl_a_debug.lib ^
  -DUDUNITS2_INCLUDE:PATH=%root%/UDUNITS-2/lib ^
  -DUDUNITS2_LIBRARY:FILE=%root%/UDUNITS-2/build/lib/Debug/udunits2.lib ^
  -DEXPAT_LIBRARY:FILE=%root%/libexpat/expat/build/Debug/expatd.lib
  msbuild nco.sln /target:build /property:configuration=debug
  if errorlevel 1 goto :eof
)

:test_nco
pushd Debug
@echo on
%root_win%\netcdf-c\build\ncgen\ncgen.exe -k netCDF-4 -b -o %root_win%\..\data\in_grp.nc %root_win%\..\data\in_grp.cdl
ncks.exe --jsn_fmt 2 -C -g g10 -v two_dmn_rec_var %root_win%\..\data\in_grp.nc
ncks.exe -v lat http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/cmap/enh/precip.mon.mean.nc
@echo off
popd
echo done