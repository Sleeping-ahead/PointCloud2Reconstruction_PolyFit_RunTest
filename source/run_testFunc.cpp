#pragma once
#include "DataStruct.h"
#include "FileProcess.h"
#include "PreProcess.h"
#include "RegionSegProcess.h"
#include "HypothesisPlaneProcess.h"
#include "ReconstructProcess.h"

namespace TempFunc
{
	//点云读取为Eigen格式
	void GetInput2Eigen(const std::string& path_input, const std::string& path_output,
		std::vector<std::vector<Eigen::Vector3f>>& cloud_input_eigen)
	{
		//获取输入文件列表
		FileProcess obj_file;
		std::vector<std::string> fileList;
		fileList = obj_file.GetFileList(path_input);

		for (const std::string& file_path : fileList)
		{
			//文件转点云
			PointCloud cloud_input;
			if (obj_file.File2PointCloud(file_path, cloud_input) == 0)
				continue;

			//转格式_eigen存储
			std::vector<Eigen::Vector3f> points(cloud_input.size());
			for (int i = 0; i < cloud_input.size(); ++i)
			{
				Eigen::Vector3f p(
					cloud_input.point(i).x(),
					cloud_input.point(i).y(),
					cloud_input.point(i).z());
				points[i] = p;
			}
			cloud_input_eigen.push_back(points);
		}

	}


	//获取目标点的k个最近邻点
	std::vector<int> GetNearPoints(const std::vector<Eigen::Vector3f>& cloud, int i_target, int k)
	{
		std::vector<std::pair<float, int>> distances;

		//计算所有点与目标点的距离
		for (int i = 0; i < cloud.size(); ++i)
		{
			//去掉目标点，避免自相关影响
			//目标点距离为零，可能会使协方差矩阵的某方向方差偏小，影响主方向估计？
			if (i != i_target)
			{
				float dist = (cloud[i_target] - cloud[i]).norm();
				distances.push_back(std::make_pair(dist, i));
			}
		}

		//按距离排序，选择最近的k个点
		std::nth_element(distances.begin(), distances.begin() + k, distances.end());
		std::vector<int> i_neighbors;
		for (int i = 0; i < k; ++i)
			i_neighbors.push_back(distances[i].second);

		return i_neighbors;
	}


	//估计法向量
	Eigen::Vector3f EstimateNormal(const std::vector<Eigen::Vector3f>& cloud, const std::vector<int>& i_neighbors)
	{
		//临近点质心
		Eigen::Vector3f p_center(0, 0, 0);
		for (int i : i_neighbors)
			p_center += cloud[i];
		p_center /= static_cast<float>(i_neighbors.size());

		//邻近点协方差矩阵
		Eigen::Matrix3f mat_covariance = Eigen::Matrix3f::Zero();
		for (int i : i_neighbors)
		{
			Eigen::Vector3f v_dis = cloud[i] - p_center;
			mat_covariance += v_dis * v_dis.transpose();
		}
		mat_covariance /= static_cast<float>(i_neighbors.size() - 1);

		//特征值分解，最小特征值对应的特征向量为法向量
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigenSlover(mat_covariance);
		Eigen::Vector3f normal = eigenSlover.eigenvectors().col(0);//对应最小特征值的特征向量

		return normal;
	}

	//保存坐标、法向量
	void SaveCloud(const std::string& path_output, const std::string& file_name,
		const std::vector<Eigen::Vector3f>& cloud,
		const std::vector<Eigen::Vector3f>& normals)
	{
		if (cloud.size() != normals.size())
		{
			std::cout << "点与法线的数目不一致" << std::endl;
			return;
		}

		std::string path_entire = path_output + file_name + ".txt";
		std::cout << "保存文件：" << path_entire << std::endl;

		std::ofstream ofs;
		ofs.flags(std::ios::fixed);
		//ofs.precision(10);//保留小数点后3位
		ofs.open(path_entire, std::ios::out);

		for (int i = 0; i < cloud.size(); ++i)
			ofs 
			<< cloud[i].x() << " "
			<< cloud[i].y() << " "
			<< cloud[i].z() << " "
			<< normals[i].x() << " "
			<< normals[i].y() << " "
			<< normals[i].z() << " "
			<< std::endl;

		ofs.close();
	}


	//输出随机颜色
	void GetRandomColor(const std::string& path_input, const std::string& path_output)
	{
		//获取输入文件列表
		FileProcess obj_file;
		std::vector<std::string> fileList;
		fileList = obj_file.GetFileList(path_input);

		for (const std::string& file_path : fileList)
		{
			//文件转点云
			PointCloud cloud_input;
			if (obj_file.File2PointCloud(file_path, cloud_input) == 0)
				continue;

			FileProcess obj_file;
			obj_file.GiveCloudRandomColor(cloud_input);
			obj_file.SaveCloud(path_output, obj_file.GetFileName(file_path) + "_color", cloud_input);

			std::cout << std::endl;
		}
	}


} //End namespace


int main_test()
{
	//设置输入、输出路径
	std::string path_input("./input/");
	std::string path_output("./output/");

	//点云计算法向量属性
	std::vector<std::vector<Eigen::Vector3f>> cloud_input_allFile;
	TempFunc::GetInput2Eigen(path_input, path_output, cloud_input_allFile);

	for (int i_file = 0; i_file < cloud_input_allFile.size(); ++i_file)
	{
		//提取点的法线信息
		std::vector<Eigen::Vector3f> cloud_input = cloud_input_allFile[i_file];
		std::vector<Eigen::Vector3f> cloud_normal(cloud_input.size());//每个点的法向量
		for (int i_p = 0; i_p < cloud_input.size(); ++i_p)
		{
			std::vector<int> i_neighbors = TempFunc::GetNearPoints(cloud_input, i_p, 15);
			Eigen::Vector3f p_normal = TempFunc::EstimateNormal(cloud_input, i_neighbors);
			cloud_normal[i_p] = p_normal;
		}

		//保存
		TempFunc::SaveCloud(path_output, "normal_" + std::to_string(i_file),
			cloud_input, cloud_normal);
	}

	return 0;
}