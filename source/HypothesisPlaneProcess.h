#pragma once
#include "DataStruct.h"
#include "GeometryCal.h"

#include <CGAL/Polygonal_surface_reconstruction/internal/point_set_with_planes.h>//带平面的输入点
#include <CGAL/Polygonal_surface_reconstruction/internal/compute_confidences.h>//计算候选面的置信度


using PNI = boost::tuple<Kernel::Point_3, Kernel::Vector_3, int>;//point_normal_int//此处需要用到boost库


class HypothesisPlaneProcess
{
public:
	
	//1、改进平面（合并点云中近似平行的区域）
	void PlaneRefine();

	//修改角度阈值_判断两区域是否平行
	inline void ChangeAngle_PlanePara(double angle);

	//2、构建包围盒的多边形网络
	void BuildMeshBbox(Polygon_mesh& mesh);

	//3、构建候选平面的多边形网络
	void BuildMeshCandidate(Polygon_mesh& bbox_mesh, Polygon_mesh& mesh_candidate);

	//4、候选平面相交细分
	void SegMeshCandidate(Polygon_mesh& mesh_candidate);

	//5.1.1、候选平面简化（提取感兴趣区域的候选平面）_激进版本
	void SimpleMeshCandidate(Polygon_mesh& mesh_candidate);

	//5.1.2、候选平面简化（提取感兴趣区域的候选平面）_保守版本
	void SimpleMeshCandidate_V2(Polygon_mesh& mesh_candidate);

	//5.2、添加额外的平面属性_到顶面的距离
	void AddMeshProperty(Polygon_mesh& mesh_candidate);

	/*main、生成候选表面*/
	void BuildHypothesisPlane(const PointCloud& cloud_region, CGAL::Surface_mesh<Kernel::Point_3>& mesh_candidate);

	//获取三个面的交点集合（三元交点组）
	inline Planes_intersections GetTriplet();

	//使多边形网络平移至原坐标
	void Coordination2Original(Polygon_mesh & model, double T_x, double T_y, double T_z);


private:

	std::vector<const Kernel::Plane_3*>		_supporting_planes; //所有平面线段和边界框面（包围盒）的支撑平面
	Planes_intersections					_triplet_intersections;
	CGAL::internal::Point_set_with_planes<Kernel> *_cloud_plane;

	//结构体_边与边的相交点
	struct EdgePos
	{
		Polygon_mesh::Edge_index edge;
		const Kernel::Point_3* pos;

		EdgePos(Polygon_mesh::Edge_index e, const Kernel::Point_3* p) : edge(e), pos(p) {}
	};

	FT _theta = static_cast<FT>(CGAL_PI * 20.0 / FT(180.0));//阈值 theta，弧度制//2024.3.8


	//1.1、计算平面段 s 上距离平面 plane_cutting 小于 dist_threshold 的点的数量
	//参数 s 指向线段的指针
	//参数 plane_cutting 指向平面的指针
	//参数 dist_threshold 距离阈值
	std::size_t calNumber_points_on_plane(const CGAL::internal::Planar_segment<Kernel>* s, const Kernel::Plane_3* plane, FT dist_threshold);

	//1.2、合并平面段 s1 和 s2
	void plane_merge(CGAL::internal::Planar_segment<Kernel> *s1, CGAL::internal::Planar_segment<Kernel> *s2);

	#pragma region 候选面集细分相关，基本照搬无修改

	//修改排序
	template <typename VT>
	inline void sort_increasing(VT& v1, VT& v2, VT& v3);

	//4.1、计算多边形mesh中的面face的交点情况，交点分类为已有点+新生成点
	void compute_intersections(const Polygon_mesh& mesh, Polygon_mesh::Face_index face, const Kernel::Plane_3* plane_cutting,
		std::vector<Polygon_mesh::Vertex_index>& existing_vts, std::vector<EdgePos>& new_vts);

	//4.2、查询半边是否存在
	bool halfedge_exists(Polygon_mesh::Vertex_index v1, Polygon_mesh::Vertex_index v2, const Polygon_mesh& mesh);

	//4.3、通过插入新点来分割边
	Polygon_mesh::Halfedge_index split_edge(Polygon_mesh& mesh, const EdgePos& ep, const Kernel::Plane_3* cutting_plane);

	//4.4、使用cuting_plane剪切面并返回新面
	std::vector<Polygon_mesh::Face_index> split_plane(Polygon_mesh::Face_index face, const Kernel::Plane_3* cutting_plane, Polygon_mesh& mesh);

	#pragma endregion

};


template <typename VT>
inline void HypothesisPlaneProcess::sort_increasing(VT& v1, VT& v2, VT& v3)
{
	VT vmin = 0;
	if (v1 < v2 && v1 < v3)
		vmin = v1;
	else if (v2 < v1 && v2 < v3)
		vmin = v2;
	else
		vmin = v3;

	VT vmid = 0;
	if ((v1 > v2 && v1 < v3) || (v1 < v2 && v1 > v3))
		vmid = v1;
	else if ((v2 > v1 && v2 < v3) || (v2 < v1 && v2 > v3))
		vmid = v2;
	else
		vmid = v3;

	VT vmax = 0;
	if (v1 > v2 && v1 > v3)
		vmax = v1;
	else if (v2 > v1 && v2 > v3)
		vmax = v2;
	else
		vmax = v3;

	v1 = vmin;
	v2 = vmid;
	v3 = vmax;
}

inline void HypothesisPlaneProcess::ChangeAngle_PlanePara(double angle)
{
	FT _theta = static_cast<FT>(CGAL_PI * angle / FT(180.0));//阈值 theta，弧度制
}

inline Planes_intersections HypothesisPlaneProcess::GetTriplet()
{
	return _triplet_intersections;
}
