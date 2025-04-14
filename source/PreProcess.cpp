#include "PreProcess.h"

void PrePorcess::Coordination2Standard(PointCloud & cloud)
{
	//1、求点云外包络立方体在xyz上的极值
	CGAL::Bbox_3 bbox = CGAL::bbox_3(cloud.points().begin(), cloud.points().end());

	//2、计算偏移量
	_xOffset = (bbox.xmax() + bbox.xmin()) / 2;
	_yOffset = (bbox.ymax() + bbox.ymin()) / 2;
	_zOffset = (bbox.zmax() + bbox.zmin()) / 2;

	//3、对每个点赋予偏移量，使点云中心移动到原点
	for (int i = 0; i < cloud.size(); ++i)
	{
		Kernel::Point_3 p_temp(cloud.point(i).x() - _xOffset,
			cloud.point(i).y() - _yOffset,
			cloud.point(i).z() - _zOffset);
		cloud.point(i) = p_temp;
	}

}

void PrePorcess::Coordination2Original(PointCloud & cloud)
{
	//对每个点赋予偏移量，使点云中心移动到原点
	for (int i = 0; i < cloud.size(); ++i)
	{
		Kernel::Point_3 p_temp(cloud.point(i).x() + _xOffset,
			cloud.point(i).y() + _yOffset,
			cloud.point(i).z() + _zOffset);
		cloud.point(i) = p_temp;
	}

}

std::vector<Kernel::Point_3> PrePorcess::CreateFootPoints(PointCloud & cloud)
{
	//1、取包围盒中部到点云底部的点云，投影至二维水平面
	CGAL::Bbox_3 bbox = CGAL::bbox_3(cloud.points().begin(), cloud.points().end());
	double bound_zmid = (bbox.zmax() + bbox.zmin()) / 2;

	std::vector<Kernel::Point_2> points_proj;
	for (size_t i = 0; i < cloud.size(); ++i)
		if (cloud.point(i).z() < bound_zmid)
			points_proj.push_back(Kernel::Point_2(cloud.point(i).x(), cloud.point(i).y()));


	//2、获取底部的外轮廓点集，以二维格网方法
	std::vector<Kernel::Point_2> points_outline_2d;
	GetOutLineByGrid(points_proj, points_outline_2d);


	//3、向外轮廓内部填充点
	SortPoints(points_outline_2d);//对构成外轮廓的多边形顶点排序
	std::vector<Kernel::Point_2> points_inside = CreatePointsInPolygon(points_outline_2d);//以一定间隔的点填充多边形区域


	//4.1、向原始点云中添加底部填充点
	for (const Kernel::Point_2& p : points_inside)
		cloud.insert(Kernel::Point_3(p.x(),
									 p.y(),
									 bbox.zmin()));


	//4.2、存储轮廓点并return
	std::vector<Kernel::Point_3> points_outline_3d(points_outline_2d.size());
	for (int i = 0; i < points_outline_2d.size(); ++i)
		points_outline_3d[i] = Kernel::Point_3(points_outline_2d[i].x(),
											   points_outline_2d[i].y(),
											   bbox.zmin());

	return points_outline_3d;
}

void PrePorcess::GetOutLineByGrid(const std::vector<Kernel::Point_2>& points, std::vector<Kernel::Point_2>& points_outline)
{
	//1、构建关于点云的二维格网，获取格网中心坐标
	std::vector<GridCenter> grid_center;

	for (int i = 0; i < points.size(); ++i)
	{
		float center_x = static_cast<float>(std::floor((points[i].x() - 0) / _size_voxel));
		float center_y = static_cast<float>(std::floor((points[i].y() - 0) / _size_voxel));
		GridCenter temp(center_x, center_y);

		//将未记录的格网索引存入grid_center
		if (std::find(grid_center.begin(), grid_center.end(), temp) == grid_center.end())
			grid_center.push_back(temp);
	}


	//2、筛选外轮廓点
	std::vector<GridCenter> grid_index_copy;
	grid_index_copy.assign(grid_center.begin(), grid_center.end());

	//为规避遮挡问题，点集绕Z轴旋转数次
	GeometryCal obj_geo;
	Eigen::Vector3f rotate_axis = Eigen::Vector3f(0, 0, 1.0f);
	for (int i_rot = 0; i_rot < 180.0f / _rotate_angle; ++i_rot)
	{
		obj_geo.Coordinate_AxisRotate(grid_center, rotate_axis, _rotate_angle);//绕z轴转_rotate_angle

		//格网索引按x值由小到大排序
		std::sort(grid_center.begin(), grid_center.end(), [](const GridCenter& a, const GridCenter& b) {
			return a.x < b.x;
		});

		//以x为检索方向，逐行查找外轮廓点
		GridCenter grid_near = grid_center[0], grid_far = grid_center[0];
		int index_near = 0, index_far = 0;
		float num_col = static_cast<float>(std::floor(grid_center[0].x));//格网列计数
		for (int i = 0; i < grid_center.size(); ++i)
		{
			if (grid_center[i].x - num_col > 1.0f)//猜测：不同尺度的格网大小应该也适用1.0f
			{
				//已开始遍历下一行数据，根据index_near及index_far对grid_center进行标记
				grid_center[index_near].isOutline = true;
				grid_center[index_far].isOutline = true;

				//更新
				num_col = static_cast<float>(std::floor(grid_center[i].x));
				grid_near = grid_center[i];
				grid_far = grid_center[i];
				index_near = i;
				index_far = i;
				continue;
			}

			//查找一行数据中的最大最小值
			if (grid_near.y < grid_center[i].y)
			{
				grid_near.y = grid_center[i].y;
				index_near = i;
			}
			if (grid_far.y > grid_center[i].y)
			{
				grid_far.y = grid_center[i].y;
				index_far = i;
			}
		}
	}

	//旋转复位
	obj_geo.Coordinate_AxisRotate(grid_center, rotate_axis, -180.0f);//反向转180度

	//存储表征外轮廓的格网中心点
	std::vector<GridCenter> grid_outline;
	for (const GridCenter& grid : grid_center)
		if (grid.isOutline)
		{
			GridCenter center_raw(grid.x * _size_voxel, grid.y * _size_voxel);//还原
			grid_outline.push_back(center_raw);
		}


	//3、计算各格网中心在原始点云中的最邻近点
	points_outline.clear();
	for (const GridCenter& p_grid : grid_outline)
	{
		//计算最近点
		double dist_min = 1e7;
		int index_near = -1;
		for (int i = 0; i < points.size(); ++i)
		{
			double dist_square = (points[i].x() - p_grid.x) * (points[i].x() - p_grid.x) +
								(points[i].y() - p_grid.y) * (points[i].y() - p_grid.y);
			if (dist_square < dist_min)
			{
				dist_min = dist_square;
				index_near = i;
			}
		}

		//存在近邻点
		if (index_near != -1)
			points_outline.push_back(points[index_near]);
	}

}

void PrePorcess::SortPoints(std::vector<Kernel::Point_2>& points)
{
	//Kernel::Point_2 center = GetCenter(points);//计算重心
	Kernel::Point_2 center = GetInner(points);//计算多边形的内点
	//冒泡排序
	size_t n = points.size();
	for (int i = 0; i < n - 1; i++)
		for (int j = 0; j < n - i - 1; j++)
			if (PointCompare(points[j], points[j + 1], center))
				swap(points[j], points[j + 1]);//前一个元素大于后一个元素，交换两个元素的位置

}

Kernel::Point_2 PrePorcess::GetCenter(const std::vector<Kernel::Point_2>& points_polygon)
{
	double x = 0, y = 0;
	size_t n = points_polygon.size();
	for (int i = 0; i < n; i++)
	{
		x += points_polygon[i].x();
		y += points_polygon[i].y();
	}
	return Kernel::Point_2(x / n, y / n);
}

Kernel::Point_2 PrePorcess::GetInner(const std::vector<Kernel::Point_2>& points_polygon)
{
	//找到最左边的顶点
	int leftmost = 0;
	for (int i = 1; i < points_polygon.size(); ++i)
		if (points_polygon[i].x() < points_polygon[leftmost].x())
			leftmost = i;

	//查询相邻的两个顶点
	size_t prev = (leftmost - 1 + points_polygon.size()) % points_polygon.size();
	size_t next = (leftmost + 1) % points_polygon.size();

	// 计算三角形的重心
	double x = (points_polygon[prev].x() + points_polygon[leftmost].x() + points_polygon[next].x()) / 3.0;
	double y = (points_polygon[prev].y() + points_polygon[leftmost].y() + points_polygon[next].y()) / 3.0;

	return Kernel::Point_2(x, y);
}

bool PrePorcess::PointCompare(const Kernel::Point_2 & a, const Kernel::Point_2 & b, const Kernel::Point_2 & center)
{
	if (a.x() == center.x() && a.y() == center.y())	return false; // 如果第一个点与重心重合，第一个点不大于第二个点
	if (b.x() == center.x() && b.y() == center.y())	return true; // 如果第二个点与重心重合，第一个点大于第二个点
	if (a.x() - center.x() >= 0 && b.x() - center.x() < 0)	return true;  // 如果第一个点在右半平面，第二个点在左半平面, 第一个点大于第二个点
	if (a.x() - center.x() < 0 && b.x() - center.x() >= 0)	return false; // 如果第一个点在左半平面，第二个点在右半平面, 第一个点不大于第二个点

	if (a.x() - center.x() == 0 && b.x() - center.x() == 0) { // 如果两个点在同一条垂直线上
		if (a.y() - center.y() >= 0 || b.y() - center.y() >= 0) { // 如果两个点都在上半平面或者都在下半平面
			return a.y() > b.y(); // 比较纵坐标大小
		}
		return b.y() > a.y(); // 否则反转纵坐标大小比较结果
	}

	// 计算两个向量的叉积
	double det = (a.x() - center.x()) * (b.y() - center.y()) - (b.x() - center.x()) * (a.y() - center.y());
	if (det < 0)
		return true; // 叉积为负，a在b的逆时针方向

	if (det > 0)
		return false; // 叉积为正，a在b的顺时针方向

	//若叉积为零，计算两个点到重心的距离
	double d1 = (a.x() - center.x()) * (a.x() - center.x()) + (a.y() - center.y()) * (a.y() - center.y());
	double d2 = (b.x() - center.x()) * (b.x() - center.x()) + (b.y() - center.y()) * (b.y() - center.y());
	return d1 > d2; // 比较距离大小
}

std::vector<Kernel::Point_2> PrePorcess::CreatePointsInPolygon(const std::vector<Kernel::Point_2>& points_polygon)
{
	std::vector<Kernel::Point_2> point_train;// 存储生成的点

	//计算多边形包围盒，赋初值
	CGAL::Bbox_2 bbox = CGAL::bbox_2(points_polygon.begin(), points_polygon.end());
	double x = bbox.xmin(), y = bbox.ymin();
	
	//生成格网中的点
	std::vector<Kernel::Point_2> points_inside;
	while (x < bbox.xmax())
	{
		while (y < bbox.ymax())
		{
			points_inside.push_back({ x, y });
			y += _size_interval;
		}
		x += _size_interval;
		y = bbox.ymin();
	}

	//判断点是否在多边形内部
	GeometryCal obj_geo;
	for (Kernel::Point_2 point : points_inside)
		if (obj_geo.isPointInPolygon(point, points_polygon))
			point_train.push_back(point);//如果点在多边形内部，将点加入到数组中

	return point_train; // 返回数组
}

