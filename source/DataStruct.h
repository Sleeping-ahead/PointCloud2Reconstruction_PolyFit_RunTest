#pragma once
#include <string>
#include <vector>

#include <Eigen/dense>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>//double内核
#include <CGAL/Point_set_3.h>//point_set类型
#include <CGAL/Polygonal_surface_reconstruction.h>//多边形曲面重建
#include <CGAL/Timer.h>//计时

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = Kernel::FT;
using PointCloud = CGAL::Point_set_3<Kernel::Point_3>;
using Polygon_mesh = CGAL::Surface_mesh<Kernel::Point_3>;
using Planes_intersections = typename std::unordered_map<const Kernel::Plane_3*,
	std::unordered_map<const  Kernel::Plane_3*,
	std::unordered_map<const Kernel::Plane_3*,
	const Kernel::Point_3*>>>;//三平面的相交点map



//格网中心点
struct GridCenter
{
	float x;//对应Col
	float y;//对应Row
	bool isOutline = false;//外轮廓标记

	GridCenter(float x, float y) : x(x), y(y) {};

	bool operator==(const GridCenter& other) const {
		return std::abs(x - other.x) <= 1e-6 && std::abs(y - other.y) <= 1e-6;
	}//操作符定义 ==

};



