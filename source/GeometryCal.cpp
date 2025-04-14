#include "GeometryCal.h"

void GeometryCal::Coordinate_AxisRotate(std::vector<Kernel::Point_3>& cloud, const Eigen::Vector3d & axis, const double & angle)
{
	double radian = angle * EIGEN_PI / 180.0f;
	for (Kernel::Point_3& p : cloud)
	{
		Eigen::AngleAxisd rotation(radian, axis);//旋转向量
		Eigen::Vector3d vec(p.x(), p.y(), p.z());
		vec = rotation * vec;
		p = Kernel::Point_3(vec.x(), vec.y(), vec.z());
	}
}

void GeometryCal::Coordinate_AxisRotate(std::vector<GridCenter>& grids, const Eigen::Vector3f& axis, const float& angle)
{
	float radian = static_cast<float>(angle * EIGEN_PI / 180.0f);
	for (GridCenter& p_center : grids)
	{
		Eigen::AngleAxisf rotation(radian, axis);//旋转向量
		Eigen::Vector3f vec(p_center.x, p_center.y, 0.0f);
		vec = rotation * vec;

		p_center.x = vec.x();
		p_center.y = vec.y();
	}
}

void GeometryCal::Coordination2Original(Polygon_mesh & model, double T_x, double T_y, double T_z)
{
	// 遍历所有的顶点并应用平移
	for (Polygon_mesh::Vertex_index v : model.vertices())
	{
		Kernel::Point_3& p = model.point(v);
		p = Kernel::Point_3(p.x() + T_x, p.y() + T_y, p.z() + T_z);
	}
}

inline void GeometryCal::GetSlopeIntercept(const Kernel::Point_2 & p1, const Kernel::Point_2 & p2, double & k, double & b)
{
	if (p1.x() == p2.x())
	{
		//如果两点相同，斜率设为无穷大，截距横坐标
		k = 1e9;
		b = p1.x();
	}
	else
	{
		//计算斜率和截距
		k = (p2.y() - p1.y()) / (p2.x() - p1.x());
		b = p1.y() - k * p1.x();
	}
}

bool GeometryCal::isLineIntersect(const Kernel::Point_2 & p1, const Kernel::Point_2 & p2, const Kernel::Point_2 & q1, const Kernel::Point_2 & q2)
{
	//计算两线段的斜率和截距
	double k1, b1, k2, b2;
	GetSlopeIntercept(p1, p2, k1, b1);
	GetSlopeIntercept(q1, q2, k2, b2);

	if (k1 == k2)
		return false;//两线段平行
	else
	{
		//计算交点横坐标，判断其是否在两条线段上
		double x = (b2 - b1) / (k1 - k2);
		return x >= std::min(p1.x(), p2.x()) && x <= std::max(p1.x(), p2.x()) &&
			x >= std::min(q1.x(), q2.x()) && x <= std::max(q1.x(), q2.x());
	}
}

bool GeometryCal::isPointInPolygon(const Kernel::Point_2 & point, const std::vector<Kernel::Point_2>& polygon)
{
	size_t n = polygon.size();//多边形的顶点个数
	int count = 0;//交点个数

	for (int i = 0; i < n; i++)
	{
		//获取边起点和终点
		Kernel::Point_2 p1 = polygon[i];
		Kernel::Point_2 p2 = polygon[(i + 1) % n];
		if (point.y() == p1.y() && point.y() == p2.y() && //点在水平边上
			point.x() >= std::min(p1.x(), p2.x()) && point.x() <= std::max(p1.x(), p2.x()))
			return true;//点在多边形内部

		if (point.y() > std::min(p1.y(), p2.y()) && point.y() <= std::max(p1.y(), p2.y()))
			if (isLineIntersect(point, Kernel::Point_2(point.x() + 1e9, point.y()), p1, p2))
				count++;//点的纵坐标在边的纵坐标范围内+水平线与边相交
	}

	return count % 2 == 1; // 如果交点个数是奇数，点在多边形内部，否则不在
}
