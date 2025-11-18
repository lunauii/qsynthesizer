cmake_minimum_required(VERSION 3.16)

project(QSynthesizer VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Qt6 REQUIRED COMPONENTS Quick)

qt_standard_project_setup(REQUIRES 6.8)

qt_add_executable(appQSynthesizer
    main.cpp
)

qt_add_qml_module(appQSynthesizer
    URI QSynthesizer
    VERSION 1.0
    QML_FILES
        RESOURCES assets/svg/cube.svg assets/svg/folder.svg assets/svg/gear.svg assets/svg/plus.svg assets/svg/saw.svg assets/svg/sine.svg assets/svg/square.svg
        RESOURCES assets/fonts/boxicons.ttf
        QML_FILES qml/colorTest.qml qml/Main.qml
        QML_FILES qml/Main.qml qml/colorTest.qml
SOURCES thirdparty/swsmodel/include/_coder_swsmodel_api.h thirdparty/swsmodel/include/_coder_swsmodel_info.h thirdparty/swsmodel/include/_coder_swsmodel_mex.h thirdparty/swsmodel/include/abs.h thirdparty/swsmodel/include/angle.h thirdparty/swsmodel/include/anynan.h thirdparty/swsmodel/include/casyi.h thirdparty/swsmodel/include/cmlri.h thirdparty/swsmodel/include/coder_array.h thirdparty/swsmodel/include/coder_array_mex.h thirdparty/swsmodel/include/fft.h thirdparty/swsmodel/include/FFTImplementationCallback.h thirdparty/swsmodel/include/filter.h thirdparty/swsmodel/include/firls.h thirdparty/swsmodel/include/gammaln.h thirdparty/swsmodel/include/hanning.h thirdparty/swsmodel/include/ifft.h thirdparty/swsmodel/include/inv.h thirdparty/swsmodel/include/kaiser.h thirdparty/swsmodel/include/length.h thirdparty/swsmodel/include/log.h thirdparty/swsmodel/include/lpcfit.h thirdparty/swsmodel/include/mean.h thirdparty/swsmodel/include/mldivide.h thirdparty/swsmodel/include/padarray.h thirdparty/swsmodel/include/qrsolve.h thirdparty/swsmodel/include/resample.h thirdparty/swsmodel/include/ResampleParser.h thirdparty/swsmodel/include/roots.h thirdparty/swsmodel/include/rt_defines.h thirdparty/swsmodel/include/rt_nonfinite.h thirdparty/swsmodel/include/rtGetInf.h thirdparty/swsmodel/include/rtGetNaN.h thirdparty/swsmodel/include/rtwtypes.h thirdparty/swsmodel/include/sort.h thirdparty/swsmodel/include/sortIdx.h thirdparty/swsmodel/include/swsmodel.h thirdparty/swsmodel/include/swsmodel_data.h thirdparty/swsmodel/include/swsmodel_initialize.h thirdparty/swsmodel/include/swsmodel_rtwutil.h thirdparty/swsmodel/include/swsmodel_terminate.h thirdparty/swsmodel/include/swsmodel_types.h thirdparty/swsmodel/include/toeplitz.h thirdparty/swsmodel/include/uniformResampleKernel.h thirdparty/swsmodel/include/upfirdnCoreImpl.h thirdparty/swsmodel/include/xcorr.h thirdparty/swsmodel/include/xdlahqr.h thirdparty/swsmodel/include/xdlanv2.h thirdparty/swsmodel/include/xgeev.h thirdparty/swsmodel/include/xnrm2.h thirdparty/swsmodel/include/xzgebal.h thirdparty/swsmodel/include/xzlarfg.h thirdparty/swsmodel/include/xzlascl.h thirdparty/swsmodel/src/_coder_swsmodel_api.cpp thirdparty/swsmodel/src/_coder_swsmodel_info.cpp thirdparty/swsmodel/src/_coder_swsmodel_mex.cpp thirdparty/swsmodel/src/abs.cpp thirdparty/swsmodel/src/angle.cpp thirdparty/swsmodel/src/anynan.cpp thirdparty/swsmodel/src/casyi.cpp thirdparty/swsmodel/src/cmlri.cpp thirdparty/swsmodel/src/fft.cpp thirdparty/swsmodel/src/FFTImplementationCallback.cpp thirdparty/swsmodel/src/filter.cpp thirdparty/swsmodel/src/firls.cpp thirdparty/swsmodel/src/gammaln.cpp thirdparty/swsmodel/src/hanning.cpp thirdparty/swsmodel/src/ifft.cpp thirdparty/swsmodel/src/inv.cpp thirdparty/swsmodel/src/kaiser.cpp thirdparty/swsmodel/src/length.cpp thirdparty/swsmodel/src/log.cpp thirdparty/swsmodel/src/lpcfit.cpp thirdparty/swsmodel/src/mean.cpp thirdparty/swsmodel/src/mldivide.cpp thirdparty/swsmodel/src/padarray.cpp thirdparty/swsmodel/src/qrsolve.cpp thirdparty/swsmodel/src/resample.cpp thirdparty/swsmodel/src/roots.cpp thirdparty/swsmodel/src/rt_nonfinite.cpp thirdparty/swsmodel/src/rtGetInf.cpp thirdparty/swsmodel/src/rtGetNaN.cpp thirdparty/swsmodel/src/sort.cpp thirdparty/swsmodel/src/sortIdx.cpp thirdparty/swsmodel/src/swsmodel.cpp thirdparty/swsmodel/src/swsmodel_data.cpp thirdparty/swsmodel/src/swsmodel_initialize.cpp thirdparty/swsmodel/src/swsmodel_rtwutil.cpp thirdparty/swsmodel/src/swsmodel_terminate.cpp thirdparty/swsmodel/src/toeplitz.cpp thirdparty/swsmodel/src/uniformResampleKernel.cpp thirdparty/swsmodel/src/upfirdnCoreImpl.cpp thirdparty/swsmodel/src/xcorr.cpp thirdparty/swsmodel/src/xdlahqr.cpp thirdparty/swsmodel/src/xdlanv2.cpp thirdparty/swsmodel/src/xgeev.cpp thirdparty/swsmodel/src/xnrm2.cpp thirdparty/swsmodel/src/xzgebal.cpp thirdparty/swsmodel/src/xzlarfg.cpp thirdparty/swsmodel/src/xzlascl.cpp
RESOURCES thirdparty/swsmodel/lib/swsmodel.a thirdparty/swsmodel/_clang-format
)

# Qt for iOS sets MACOSX_BUNDLE_GUI_IDENTIFIER automatically since Qt 6.1.
# If you are developing for iOS or macOS you should consider setting an
# explicit, fixed bundle identifier manually though.
set_target_properties(appQSynthesizer PROPERTIES
#    MACOSX_BUNDLE_GUI_IDENTIFIER com.example.appQSynthesizer
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    MACOSX_BUNDLE TRUE
    WIN32_EXECUTABLE TRUE
)

target_link_libraries(appQSynthesizer
    PRIVATE Qt6::Quick
)

include(GNUInstallDirs)
install(TARGETS appQSynthesizer
    BUNDLE DESTINATION .
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
