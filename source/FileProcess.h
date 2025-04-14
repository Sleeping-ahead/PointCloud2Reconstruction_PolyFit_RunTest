#pragma once
#include "DataStruct.h"

#include <iostream>
#include <fstream>
#include <sstream>
#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING
#include <experimental/filesystem>//C++17 文件操作
//获取系统时间，仅用于打印输出时间
//#pragma warning(disable: 4996)//ctime unsafe
//#include <chrono>
//#include <ctime>

#include <CGAL/IO/read_points.h>//读取
//#include <CGAL/IO/read_ply_points.h>//读取ply文件


class FileProcess
{
public:

	//获取文件名列表
	std::vector<std::string>GetFileList(const std::string& path);

	//返回文件名
	std::string GetFileName(const std::string& path_file);

	//加载文件为点云
	int File2PointCloud(const std::string& path_input, PointCloud& cloud);

	//赋予（单个）点云随机颜色
	void GiveCloudRandomColor(PointCloud& cloud);

	//保存点云文件_Kernel::Point_3
	void SaveCloud(const std::string& path_output, const std::string& file_name,
		const std::vector<Kernel::Point_3>& data);

	//保存点云文件_CGAL::Point_set_3<Kernel::Point_3>
	void SaveCloud(const std::string& path_output, const std::string& file_name,
		const PointCloud& data);

	//输出无色的obj文件
	bool SaveOBJ(const std::string & path_output, const std::string & file_name,
		const CGAL::Surface_mesh<Kernel::Point_3>& mesh);

	//输出随机颜色的obj文件
	bool SaveColorOBJ(const std::string& path_output, const std::string& file_name,
		const CGAL::Surface_mesh<Kernel::Point_3>& mesh);

	//输出随机颜色的off文件//可弃用，格式的显示比较单一
	bool SaveColorOFF(const std::string& path_output, const std::string& file_name,
		const CGAL::Surface_mesh<Kernel::Point_3>& mesh);


private:


};

