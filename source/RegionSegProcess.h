#pragma once
#include "DataStruct.h"

#include <CGAL/pca_estimate_normals.h>//pca计算法向量
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>//区域增长算法
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>


using PointVectorPair = std::pair<Kernel::Point_3, Kernel::Vector_3>;

using Point_map = PointCloud::Point_map;
using Normal_map = PointCloud::Vector_map;
using PNI = boost::tuple<Kernel::Point_3, Kernel::Vector_3, int>;//point_normal_int
using Neighbor_query = CGAL::Shape_detection::Point_set::K_neighbor_query<Kernel, PointCloud, Point_map>;//k近邻查询
using Region_type = CGAL::Shape_detection::Point_set::Least_squares_plane_fit_region<Kernel, PointCloud, Point_map, Normal_map>;
using Region_growing = CGAL::Shape_detection::Region_growing<PointCloud, Neighbor_query, Region_type>;


class RegionSegProcess
{
public:

	//PCA法估算点云法向量
	void EstimateNormals(const PointCloud& cloud, PointCloud& cloud_normal);

	//基于区域增长的点云分割，rueturn:未满足条件的区域
	PointCloud RegionGrowing(const PointCloud& cloud_normal, PointCloud& cloud_seg);

	//修改比例参数
	inline void ChangeParameter(const int& num_kSearch = 18, const FT& max_distance_point2plane = 2,
		const FT& max_accepted_angle = 30, const int& min_region_pointNum = 250);

	//查询较大跨幅点云//temp
	void FindBigVolumeCloud(const PointCloud& cloud_seg, std::vector<PointCloud>& cloud_bigVolume);


private:

	int _num_kSearch = 18; //k近邻搜索的邻域点数

	//稀疏点云小建筑的参数
	FT	_max_distance_point2plane = FT(2);//距离（拟合）平面的最大距离
	FT	_max_accepted_angle = FT(30);//法线相似性的最大角度
	int	_min_region_pointNum = 250;//集群的最小值
	//密集点云大建筑的参数
	//FT	_max_distance_point2plane = FT(2);//距离（拟合）平面的最大距离
	//FT	_max_accepted_angle = FT(30);//法线相似性的最大角度
	//int	_min_region_pointNum = 1000;//集群的最小值

};

inline void RegionSegProcess::ChangeParameter(const int& num_kSearch, const FT& max_distance_point2plane,
	const FT& max_accepted_angle, const int& min_region_pointNum)
{
	_num_kSearch = num_kSearch;
	_max_distance_point2plane = max_distance_point2plane;
	_max_accepted_angle = max_accepted_angle;
	_min_region_pointNum = min_region_pointNum;
}
