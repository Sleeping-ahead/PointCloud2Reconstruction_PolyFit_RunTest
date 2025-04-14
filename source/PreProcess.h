#pragma once
#include "DataStruct.h"
#include "GeometryCal.h"

//#include <CGAL/intersections.h>//求交计算

//#include <CGAL/Alpha_shape_2.h>//alpha凹包计算
//#include <CGAL/Alpha_shape_vertex_base_2.h>
//#include <CGAL/Alpha_shape_face_base_2.h>
//#include <CGAL/Delaunay_triangulation_2.h>


//Alpha形状提取
//using Vb = CGAL::Alpha_shape_vertex_base_2<Kernel>;// Alpha形状的顶点基类类型
//using Fb = CGAL::Alpha_shape_face_base_2<Kernel>;// Alpha形状的面基类类型
//using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;// 三角剖分的数据结构类型
//using Triangulation_2 = CGAL::Delaunay_triangulation_2<Kernel, Tds>;// Delaunay三角剖分类型
//using Alpha_shape_2 = CGAL::Alpha_shape_2<Triangulation_2>;// Alpha形状类型



class PrePorcess
{
public:

	//原始坐标系到标准坐标系的三轴偏移量
	double _xOffset = 0.0f;
	double _yOffset = 0.0f;
	double _zOffset = 0.0f;


	//使点云平移至标准坐标
	void Coordination2Standard(PointCloud& cloud);

	//使点云平移至原坐标
	void Coordination2Original(PointCloud& cloud);

	//创建点云的底部填充点云，return:底部的外轮廓点集
	std::vector<Kernel::Point_3> CreateFootPoints(PointCloud& cloud);

	//参数修改
	inline void ChangeParameter(const float& size_voxel = 1.0f,
		const float& rotate_angle = 18.0f,
		const float& size_interval = 0.2f);

private:
	
	float _size_voxel = 1.0f;//二维格网大小
	float _rotate_angle = 18.0f;//每次旋转的角度大小
	float _size_interval = 0.2f;//底部填充点的间隔大小

	//通过格网获取二维点云的外轮廓点集
	void GetOutLineByGrid(const std::vector<Kernel::Point_2>& points, std::vector<Kernel::Point_2>& points_outline);

	#pragma region 多边形顶点排序

	//多边形顶点按逆时针方向排序//适用于凸多边形
	void SortPoints(std::vector<Kernel::Point_2>& points);

	//计算多边形的重心
	Kernel::Point_2 GetCenter(const std::vector<Kernel::Point_2>& points_polygon);

	//计算多边形内点
	Kernel::Point_2 GetInner(const std::vector<Kernel::Point_2>& points_polygon);

	//比较两个点的大小关系
	bool PointCompare(const Kernel::Point_2& a, const Kernel::Point_2& b, const Kernel::Point_2& center);

	#pragma endregion

	//多边形内部按间隔生成点
	std::vector<Kernel::Point_2> CreatePointsInPolygon(const std::vector<Kernel::Point_2>& points_polygon);



};

inline void PrePorcess::ChangeParameter(const float & size_voxel, const float & rotate_angle, const float & size_interval)
{
	_size_voxel = size_interval;
	_rotate_angle = rotate_angle;
	_size_interval = size_interval;
}
