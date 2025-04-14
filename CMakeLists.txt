cmake_minimum_required(VERSION 3.15)
project(cloud2Mesh)

### 设置配置类型
set(CMAKE_CONFIGURATION_TYPES "Release" CACHE STRING "Limited configurations" FORCE)


### 指定外部依赖库的查找目录
# 注：只需修改每个三方库的根目录
set(BOOST_DIR_INCLUDE "D:/BOOST/include/boost_1_82_0" CACHE PATH "Boost include directory")
set(CGAL_DIR_ROOT "D:/CGAL/CGAL-5.5.2" CACHE PATH "CGAL root directory")
set(SCIP_DIR_ROOT "D:/SCIP/SCIPOptSuite 8.0.3" CACHE PATH "SCIP root directory")
set(EIGEN_DIR_INCLUDE "D:/eigen-3.4.0" CACHE PATH "Eigen include directory")

# 检测路径有效性
function(check_path_exist lib_name lib_path)
    if(NOT EXISTS ${lib_path})
        message(FATAL_ERROR "The specified ${lib_name} does not exist: ${lib_path}")
    endif()
endfunction()

check_path_exist("BOOST_DIR_INCLUDE" "${BOOST_DIR_INCLUDE}")
check_path_exist("CGAL_DIR_ROOT" "${CGAL_DIR_ROOT}")
check_path_exist("SCIP_DIR_ROOT" "${SCIP_DIR_ROOT}")
check_path_exist("EIGEN_DIR_INCLUDE" "${EIGEN_DIR_INCLUDE}")

# 查找外部库依赖项
# CGAL
set(CGAL_DIR_INCLUDE 
    "${CGAL_DIR_ROOT}/include"
    "${CGAL_DIR_ROOT}/auxiliary/gmp/include"
)
set(CGAL_DIR_LIB "${CGAL_DIR_ROOT}/auxiliary/gmp/lib")
set(CGAL_LIB
    "${CGAL_DIR_ROOT}/auxiliary/gmp/lib/libgmp-10.lib"
    "${CGAL_DIR_ROOT}/auxiliary/gmp/lib/libmpfr-4.lib"
)
# SCIP
set(SCIP_DIR_INCLUDE "${SCIP_DIR_ROOT}/include")
set(SCIP_DIR_LIB "${SCIP_DIR_ROOT}/lib")
set(SCIP_LIB "${SCIP_DIR_LIB}/libscip.lib")

# 查找自定义文件
file(GLOB _CPP ${PROJECT_SOURCE_DIR}/source/*.cpp) # 搜索源文件
file(GLOB _H ${PROJECT_SOURCE_DIR}/source/*.h) # 搜索头文件

# 指定包含文件路径
include_directories(
    ${CGAL_DIR_INCLUDE}
    ${SCIP_DIR_INCLUDE}
    ${EIGEN_DIR_INCLUDE}
    ${BOOST_DIR_INCLUDE}
) 
# 指定静态库路径
link_directories(
    ${CGAL_DIR_LIB}
    ${SCIP_DIR_LIB}
) 

# dll文件拷贝
# set(DLL_SOURCE_DIR ${PROJECT_SOURCE_DIR}/dll) # DLL文件源路径
set(DLL_DESTINATION_DIR ${PROJECT_BINARY_DIR}/Release) # DLL文件目标路径
file(GLOB DLL_SOURCE_DIR ${PROJECT_SOURCE_DIR}/dll/*.dll) # 搜索DLL文件
file(MAKE_DIRECTORY ${DLL_DESTINATION_DIR})
file(COPY ${DLL_SOURCE_DIR} # 指定复制源路径
     DESTINATION ${DLL_DESTINATION_DIR} # 指定目标路径
     FILES_MATCHING PATTERN "*.dll") # FILES_MATCHING只复制符合特定模式的文件；PATTERN定义匹配的后缀名


### 构建项目
# 创建输入、输出文件夹
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/input) # CMAKE_BINARY_DIR指向执行cmake命令的目录
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)

# 定义可执行文件
set(_PROJECT_LIST ${_H} ${_CPP})
add_executable(${PROJECT_NAME} ${_PROJECT_LIST}) # 构建程序,记得在_PROJECT_LIST包含.h，否则其变为外部依赖项
target_link_libraries(${PROJECT_NAME}
${CGAL_LIB}
${SCIP_LIB}
) #链接库文件


### 组织显示逻辑
# CMake项目放入文件夹CMakePredefinedTargets
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "CMakePredefinedTargets")

# 指定源文件和头文件的逻辑文件夹
foreach(FILE ${_PROJECT_LIST}) 
    # 获取文件目录
    get_filename_component(PARENT_DIR "${FILE}" DIRECTORY)

    # 文件分组为 "Source Files"和"Header Files"
    if ("${FILE}" MATCHES ".*\\.cpp")
       set(GROUP "Source Files")
    elseif("${FILE}" MATCHES ".*\\.h")
       set(GROUP "Header Files")
    endif()

    source_group("${GROUP}" FILES "${FILE}") # 指定组名，将对象添加到分组
endforeach()

