/*
@inproceedings{nan2017polyfit,
  title={Polyfit: Polygonal surface reconstruction from point clouds},
  author={Nan, Liangliang and Wonka, Peter},
  booktitle={Proceedings of the IEEE International Conference on Computer Vision},
  pages={2353--2361},
  year={2017}
}
*/

#pragma once
#include "DataStruct.h"
#include "FileProcess.h"
#include "PreProcess.h"
#include "RegionSegProcess.h"
#include "HypothesisPlaneProcess.h"
#include "ReconstructProcess.h"

//点云到多边形网络的工作流
void PointCloud2Mesh(const std::string& path_input, const std::string& path_output)
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
	
#if 1
		PrePorcess obj_pre;
		std::vector<Kernel::Point_3> point_outline;
		obj_pre.Coordination2Standard(cloud_input);//点云平移至坐标原点
		//obj_pre.ChangeParameter();
		point_outline = obj_pre.CreateFootPoints(cloud_input);//点云底部填充		
		//obj_file.SaveCloud(path_output, obj_file.GetFileName(file_path) + "_outline", point_outline);
		//obj_file.SaveCloud(path_output, obj_file.GetFileName(file_path) + "_foot", cloud_input);

#endif // 点云预处理

#if 1
		RegionSegProcess obj_seg;
		PointCloud cloud_normal;
		//obj_seg.ChangeParameter(18, 2.5, 45, 350);//参数修改
		//obj_seg.ChangeParameter(12, 2, 30, 50);//参数修改
		obj_seg.EstimateNormals(cloud_input, cloud_normal);//估计点云的法向量

		PointCloud cloud_seg;
		obj_seg.RegionGrowing(cloud_normal, cloud_seg);//区域生长分割
		//弃用，原数据平移后可能因位数问题导致法向量偏移，致使重建面片偏移
		//obj_pre.Coordination2Original(cloud_seg);//分割完毕，将数据平移至原始坐标
		PointCloud cloud_seg_copy(cloud_seg);
		obj_pre.Coordination2Original(cloud_seg_copy);//分割完毕，将数据平移至原始坐标
		obj_file.SaveCloud(path_output, obj_file.GetFileName(file_path) + "_region", cloud_seg);

#endif // 法向量计算、区域分割

#if 1
		HypothesisPlaneProcess obj_hp;
		CGAL::Surface_mesh<Kernel::Point_3> mesh_candidate;//多边形网络的候选表面集
		obj_hp.BuildHypothesisPlane(cloud_seg, mesh_candidate);//构建点云的面片集合
		//obj_hp.Coordination2Original(mesh_candidate, obj_pre._xOffset, obj_pre._yOffset, obj_pre._zOffset);//网络平移至原始坐标
		obj_file.SaveColorOBJ(path_output, obj_file.GetFileName(file_path) + "_mesh_candidate", mesh_candidate);
		//obj_file.SaveOBJ(path_output, obj_file.GetFileName(file_path) + "_mesh_candidate_nocolor", mesh_candidate);

		ReconstructProcess obj_recon;
		CGAL::Surface_mesh<Kernel::Point_3> mesh_model;//模型
		//obj_recon.UseExtraFactor(0.05f);
		//obj_recon.Reconstruct(mesh_candidate, mesh_model, obj_hp.GetTriplet(), 0.43f, 0.27f, 0.25f);//筛选构成模型的面片集合
		obj_recon.Reconstruct(mesh_candidate, mesh_model, obj_hp.GetTriplet(), 0.43f, 0.27f, 0.30f);//默认的推荐参数
		obj_file.SaveOBJ(path_output, obj_file.GetFileName(file_path) + "_model", mesh_model);
		obj_file.SaveColorOBJ(path_output, obj_file.GetFileName(file_path) + "_model_color", mesh_model);


		obj_recon.MergePlane(mesh_model);//合并细碎面片，简化模型面数
		obj_file.SaveColorOBJ(path_output, obj_file.GetFileName(file_path) + "_model_simple", mesh_model);

#endif // 生成候选平面、重建模型

		std::cout << std::endl;
	}
}


//int main()
int main()
{
	//设置输入、输出路径
	std::string path_input("./input/");
	std::string path_output("./output/");

	//点云生成模型
	PointCloud2Mesh(path_input, path_output);

	return 0;
}