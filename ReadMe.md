### 项目描述
照搬NanLiangLiang副教授于2017年提出的PolyFit方法，在此基础上做了点改进尝试和代码注释，整个程序可分为下述四步骤：   
（1）预处理：几何中心归零；点云底面填充。   
（2）结合平面检测的区域增长聚类。   
（3）区域点云的平面拟合；平面相交细分；剪除冗余面。   
（4）建立混合整数规划最小化问题，以边、面为变量构造能量函数；用求解器筛选出最优的面片组合；面片整合。

### PolyFit引用声明
>@inproceedings{nan2017polyfit,   
  title={Polyfit: Polygonal surface reconstruction from point clouds},   
  author={Nan, Liangliang and Wonka, Peter},   
  booktitle={Proceedings of the IEEE International Conference on Computer Vision},   
  pages={2353--2361},   
  year={2017}   
}

### 环境与平台
C++语言编写，windows环境，编译平台VS2017。依赖项：BOOST 1_82、[CGAL 5.5.2](https://github.com/CGAL/cgal/releases/download/v5.5.2/CGAL-5.5.2.zip)、[Eigen 3.4.0](https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip)、SCIP 8.0.3。   
这是当初可编译的依赖环境，更新的版本应该也能兼容。像[BOOST 1_85](https://boostorg.jfrog.io/artifactory/main/release/1.85.0/source/boost_1_85_0.7z)、[SCIP 9.0.0](https://scipopt.org/download.php?fname=SCIPOptSuite-9.0.0-win64-VS15.exe)。

### 构建项目
1、下载好程序的源码，解压，并创建一个构建项目的文件夹“build”。

2、打开解压路径下的CmakeLists.txt，将“指定外部依赖库的查找目录”下面的四个库目录修改为自己电脑上的。

3、两种构建方法，任选一种即可：   
（1）Cmake程序构建。打开cmake，设置源代码路径为解压路径，构建文件路径为刚创建的build文件夹。点击【 configure 】，设置相应的开发环境，设置生成平台为x64，点击【Finish】后静待构建完毕。

（2）命令行构建。打开cmd命令行，先检查下是否把cmake加到环境变量，然后cd到build文件夹，输入
>cmake -G "Visual Studio 15 2017" -A x64 ..

静待构建完毕。

### 运行
1、输入数据：点入build文件夹下的input文件夹，放入点云的三维空间坐标数据，.txt、.xyz或.ply格式都行，可以批处理，data_test文件夹内有示例数据。

2、运行：双击build文件夹下的.sln文件，待编译器解析就绪后，在解决方案资源管理器中右键项目“cloud2Mesh”，将其设为启动项目，运行。

正常情况下等命令行跑完就行，输出结果默认在input文件夹下。

可以右键项目->属性，查看附加包含目录、附加库目录和附加依赖项存不存在，如果没有的话可能是前面CmakeLists.txt的四个库路径指定错了，可以自己手动把缺的补上。

3、输出结果：程序跑完，每一份输入会得到4份输出，位于output文件夹下：与点云数据同名的聚类点云.ply、细分的候选平面集合.obj、未赋色的建筑物模型和随机颜色的建筑物模型。

### 更多说明
可访问我的B站空间查找相关视频，精力有限仅录制4小节。   
https://space.bilibili.com/179235151   
或访问挂载在github上的博客，翻阅相应的文章，包含整个流程的梳理及代码注解。   
https://sleeping-ahead.github.io/categories/点云三维重建/

项目并不完美（this is not a perfect job），其中的谬误、不足之处，文中也有所提及。水平有限，大佬们海涵。   
希望后继者更快上手，少走些弯路XD