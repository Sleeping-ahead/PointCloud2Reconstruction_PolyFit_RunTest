#include "RegionSegProcess.h"

void RegionSegProcess::EstimateNormals(const PointCloud & cloud, PointCloud & cloud_normal)
{
	//1、转格式，cloud -> list_normal
	std::list<PointVectorPair> list_normal;//二元列表_法向量
	for (int i = 0; i < cloud.size(); ++i)
	{
		PointVectorPair pair_temp;
		pair_temp.first = cloud.point(i);
		list_normal.push_back(pair_temp);
	}

	
	//2、法向量计算
	//使用固定邻近点来计算法线，pca_estimate_normals()需要一个范围的点，以及属性映射来访问每个点的位置和法线
	CGAL::pca_estimate_normals<CGAL::Parallel_if_available_tag>(list_normal,
		_num_kSearch,
		CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())
		.normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));

	///输出测试_list格式的法向量
	//for (const PointVectorPair& pair : list_normal)
	//	std::cout << pair.first << "  " << pair.second << std::endl;


#if 0
	//3、可选操作：法线定向。//后续用不上
	//需要包含#include <CGAL/mst_orient_normals.h>
	//mst_orient_normals()需要一系列点，以及属性映射来访问每个点的位置和法线。
	std::list<PointVectorPair>::iterator unoriented_points_begin =
		CGAL::mst_orient_normals(list_normal, _num_kSearch,
			CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())
			.normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));
	//删除具有无方向法线的点，计划调用期望定向法线的重建算法
	list_normal.erase(unoriented_points_begin, list_normal.end());
#endif //3、法线定向，后续用不上


	//4、转格式，list_normal -> cloud_normal
	cloud_normal.clear();
	cloud_normal.add_normal_map();
	for (const PointVectorPair& pair : list_normal)
		cloud_normal.insert(pair.first, pair.second);

}

PointCloud RegionSegProcess::RegionGrowing(const PointCloud & cloud_normal, PointCloud & cloud_seg)
{
	//1、创建类的实例
	Neighbor_query neighbor_query(
		cloud_normal,
		_num_kSearch,
		cloud_normal.point_map());//k近邻搜索邻域查询的实例
	//Neighbor_query neighbor_query(
	//	cloud_normal,
	//	search_sphere_radius);//半径搜索邻域查询的实例//易报错，弃用
	Region_type region_type(
		cloud_normal,
		_max_distance_point2plane, _max_accepted_angle, _min_region_pointNum,
		cloud_normal.point_map(), cloud_normal.normal_map());//区域类型的实例
	Region_growing region_growing(
		cloud_normal,
		neighbor_query, region_type);//区域增长的实例


	//2、区域增长分割算法
	CGAL::Timer timer;//时间计数的实例
	timer.start();
	std::vector<std::vector<std::size_t>> regions;
	region_growing.detect(std::back_inserter(regions));//算法运行
	timer.stop();
	std::cout << "有 " << regions.size() << " 个分割区域已生成，运行耗时 " << timer.time() << " 秒" << std::endl;


	//3、根据region中的区域索引，拷贝cloud_normal至cloud_seg中，附上随机颜色
	cloud_seg.clear();
	cloud_seg.add_normal_map();//为cloud_seg添加normal属性
	cloud_seg.add_property_map<int>("region_map").first;//为cloud_seg添加区域索引属性
	cloud_seg.add_property_map<unsigned char>("red").first;//为cloud_seg添加color属性
	cloud_seg.add_property_map<unsigned char>("green").first;
	cloud_seg.add_property_map<unsigned char>("blue").first;

	for (int i = 0; i < regions.size(); ++i)
	{
		//设置随机颜色
		const unsigned char r = static_cast<unsigned char>(std::rand() % 256);
		const unsigned char g = static_cast<unsigned char>(std::rand() % 256);
		const unsigned char b = static_cast<unsigned char>(std::rand() % 256);

		for (const auto& index : regions[i])
		{
			PointCloud cloud_piece;
			cloud_piece.add_normal_map();
			const auto& key = *(cloud_normal.begin() + index);//根据regions[i]建立对cloud_normal的索引
			cloud_piece.insert(cloud_normal, key);
			cloud_piece.add_property_map<int>("region_map", i).first;
			cloud_piece.add_property_map<unsigned char>("red", r).first;
			cloud_piece.add_property_map<unsigned char>("green", g).first;
			cloud_piece.add_property_map<unsigned char>("blue", b).first;

			cloud_seg.insert(cloud_piece, *cloud_piece.begin());
		}
	}


	//4、统计未分割的区域
	std::vector<std::size_t> unassigned_items;//因未满足条件而未被分割点的区域
	region_growing.unassigned_items(std::back_inserter(unassigned_items));

	PointCloud cloud_unassigned;
	for (const auto index : unassigned_items)
	{
		const auto& key = *(cloud_normal.begin() + index);//索引
		const Kernel::Point_3& point = get(cloud_normal.point_map(), key);
		cloud_unassigned.insert(point);
	}

	return cloud_unassigned;
}

void RegionSegProcess::FindBigVolumeCloud(const PointCloud& cloud_seg, std::vector<PointCloud>& cloud_bigVolume)
{
	//1、确定分割区域数量
	std::vector<int> num_region;
	std::vector<PointCloud> cloud_temp;
	for (int i = 0; i < cloud_seg.size(); ++i)
	{
		int region_index = cloud_seg.property_map<int>("region_map").first[i];
		if (std::find(num_region.begin(), num_region.end(), region_index) == num_region.end())
			num_region.push_back(region_index);
	}
	cloud_temp.resize(num_region.size());

	//2、存储点
	for (int i = 0; i < cloud_seg.size(); ++i)
	{
		int region_index = cloud_seg.property_map<int>("region_map").first[i];
		cloud_temp[region_index].insert(cloud_seg.point(i));
	}

	//3、根据包围盒体积筛选区域
	cloud_bigVolume.clear();
	for (int i = 0; i < cloud_temp.size(); ++i)
	{
		//计算包围盒
		CGAL::Bbox_3 bbox = CGAL::bbox_3(cloud_temp[i].points().begin(), cloud_temp[i].points().end());

		FT dx = bbox.xmax() - bbox.xmin();
		FT dy = bbox.ymax() - bbox.ymin();
		FT dz = bbox.zmax() - bbox.zmin();
		FT volume = dx * dy * dz;

		//筛选出较大体积的区域点云
		if (volume > 100.0f)
			cloud_bigVolume.push_back(cloud_temp[i]);

	}




}