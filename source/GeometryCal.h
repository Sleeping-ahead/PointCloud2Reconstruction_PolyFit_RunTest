#pragma once
#include "DataStruct.h"


class GeometryCal
{
public:

	//坐标变换_绕轴旋转一定角度_三维点云
	void Coordinate_AxisRotate(std::vector<Kernel::Point_3>& cloud, const Eigen::Vector3d& axis, const double& angle);

	//坐标变换_绕轴旋转一定角度_格网
	void Coordinate_AxisRotate(std::vector<GridCenter>& grids, const Eigen::Vector3f& axis, const float& angle);

	//使多边形网络平移至原坐标
	void Coordination2Original(Polygon_mesh & model, double T_x, double T_y, double T_z);


	//————————二维操作————————//

	//计算两个点的斜率和截距
	inline void GetSlopeIntercept(const Kernel::Point_2& p1, const Kernel::Point_2& p2, double& k, double& b);

	//判断两条线段是否相交
	bool isLineIntersect(const Kernel::Point_2& p1, const Kernel::Point_2& p2,
		const Kernel::Point_2& q1, const Kernel::Point_2& q2);

	//判断一个点是否在多边形内部
	bool isPointInPolygon(const Kernel::Point_2& point, const std::vector<Kernel::Point_2>& polygon);


private:

};

