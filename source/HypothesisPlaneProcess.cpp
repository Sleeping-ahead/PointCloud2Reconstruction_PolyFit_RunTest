#include "HypothesisPlaneProcess.h"

void HypothesisPlaneProcess::PlaneRefine()
{
	std::vector< CGAL::internal::Planar_segment<Kernel>* >& segments = _cloud_plane->planar_segments();//获取 point_set_ 对象的 planar_segments() 成员变量的引用

	//1、遍历所有平面段，计算每个平面段的平均最大距离
	FT avg_max_dist = 0;
	for (std::size_t i = 0; i < segments.size(); ++i)
	{
		CGAL::internal::Planar_segment<Kernel>* s = segments[i];
		const Kernel::Plane_3* plane = s->fit_supporting_plane();//拟合平面(用户可能提供无效的平面拟合

		FT max_dist = -(std::numeric_limits<FT>::max)();//初始化最大距离为无穷小
		for (std::size_t j = 0; j < s->size(); ++j)
		{
			std::size_t idx = s->at(j);//点索引
			const Kernel::Point_3& p = _cloud_plane->point_map()[idx];
			FT sdist = CGAL::squared_distance(*plane, p);//计算点到平面的距离
			max_dist = (std::max)(max_dist, std::sqrt(sdist));//更新最大距离
		}

		avg_max_dist += max_dist;//累加最大距离
	}
	//计算平均最大距离
	avg_max_dist /= segments.size();
	avg_max_dist /= FT(2.0);


	//2、平面段合并
	//FT theta = static_cast<FT>(CGAL_PI * 10.0 / FT(180.0));//阈值 theta，弧度制//放公有变量处
	bool merged = false;//是否有平面段被合并
	do {
		merged = false;

		//将平面段按照大小进行排序
		//点较少的分段可信度较低，因此应首先合并
		std::sort(segments.begin(), segments.end(), CGAL::internal::SegmentSizeIncreasing<CGAL::internal::Planar_segment<Kernel>>());

		//遍历平面段，plane1为当前，plane2为下一个
		for (std::size_t i = 0; i < segments.size(); ++i)
		{
			CGAL::internal::Planar_segment<Kernel>* s1 = segments[i];
			const Kernel::Plane_3* plane1 = s1->supporting_plane();
			Kernel::Vector_3 n1 = plane1->orthogonal_vector();
			CGAL::internal::normalize<FT, Kernel::Vector_3>(n1);

			FT num_threshold = s1->size() / FT(5.0);//阈值 num_threshold
			for (std::size_t j = i + 1; j < segments.size(); ++j)
			{
				CGAL::internal::Planar_segment<Kernel>* s2 = segments[j];
				const Kernel::Plane_3* plane2 = s2->supporting_plane();
				Kernel::Vector_3 n2 = plane2->orthogonal_vector();
				CGAL::internal::normalize<FT, Kernel::Vector_3>(n2);

				//如果两个平面段的正交向量之间的内积大于 theta，
				if (std::abs(n1 * n2) > std::cos(_theta))
				{
					//且在另一平面上的点数大于 num_threshold，合并两平面
					std::size_t set1on2 = calNumber_points_on_plane(s1, plane2, avg_max_dist);
					std::size_t set2on1 = calNumber_points_on_plane(s2, plane1, avg_max_dist);
					if (set1on2 > num_threshold || set2on1 > num_threshold)
					{
						plane_merge(s1, s2);
						merged = true;
						break;
					}
				}
			}
			if (merged)
				break;
		}
	} while (merged);//merged 若为false，则此时所有平面都不适合被合并，跳出循环

	std::sort(segments.begin(), segments.end(), CGAL::internal::SegmentSizeDecreasing<CGAL::internal::Planar_segment<Kernel>>());//平面段按大小逆序排序

	//3、存储所有的支撑平面->_supporting_planes
	for (std::size_t i = 0; i < segments.size(); ++i)
	{
		CGAL::internal::Planar_segment<Kernel>* s = segments[i];
		const Kernel::Plane_3* plane = s->supporting_plane();
		_supporting_planes.push_back(plane);
	}


}

void HypothesisPlaneProcess::BuildMeshBbox(Polygon_mesh & mesh)
{
	//1、计算包围盒半径、偏移量
	const Kernel::Iso_cuboid_3& box = CGAL::bounding_box(_cloud_plane->point_map().begin(), _cloud_plane->point_map().end());//计算点云的包围盒

	FT dx = box.xmax() - box.xmin();
	FT dy = box.ymax() - box.ymin();
	FT dz = box.zmax() - box.zmin();
	FT radius = FT(0.5) * std::sqrt(dx * dx + dy * dy + dz * dz);
	FT offset = radius * FT(0.05);


	//2、创建一个更大的包围盒，以确保所有点都包含其中
	FT xmin = box.xmin() - offset, xmax = box.xmax() + offset;
	FT ymin = box.ymin() - offset, ymax = box.ymax() + offset;
	FT zmin = box.zmin() - offset, zmax = box.zmax() + offset;


	//3、更新 mesh 中的多边形网络为点集 points 的包围盒
	mesh.clear();//清空多边形网络
	//添加包围盒的8个顶点
	Polygon_mesh::Vertex_index v0 = mesh.add_vertex(Kernel::Point_3(xmin, ymin, zmin));  // 0
	Polygon_mesh::Vertex_index v1 = mesh.add_vertex(Kernel::Point_3(xmax, ymin, zmin));  // 1
	Polygon_mesh::Vertex_index v2 = mesh.add_vertex(Kernel::Point_3(xmax, ymin, zmax));  // 2
	Polygon_mesh::Vertex_index v3 = mesh.add_vertex(Kernel::Point_3(xmin, ymin, zmax));  // 3
	Polygon_mesh::Vertex_index v4 = mesh.add_vertex(Kernel::Point_3(xmax, ymax, zmax));  // 4
	Polygon_mesh::Vertex_index v5 = mesh.add_vertex(Kernel::Point_3(xmax, ymax, zmin));  // 5
	Polygon_mesh::Vertex_index v6 = mesh.add_vertex(Kernel::Point_3(xmin, ymax, zmin));  // 6
	Polygon_mesh::Vertex_index v7 = mesh.add_vertex(Kernel::Point_3(xmin, ymax, zmax));  // 7

	//添加包围盒的6个面
	mesh.add_face(v0, v1, v2, v3);
	mesh.add_face(v1, v5, v4, v2);
	mesh.add_face(v1, v0, v6, v5);
	mesh.add_face(v4, v5, v6, v7);
	mesh.add_face(v0, v3, v7, v6);
	mesh.add_face(v2, v4, v7, v3);


	//添加多边形属性
	typename Polygon_mesh::template Property_map<Polygon_mesh::Face_index, const Kernel::Plane_3*> face_supporting_planes =
		mesh.template add_property_map<Polygon_mesh::Face_index, const Kernel::Plane_3*>("f:supp_plane").first;//每个面的支撑平面

	typename Polygon_mesh::template Property_map<Polygon_mesh::Edge_index, std::set<const Kernel::Plane_3*> > edge_supporting_planes
		= mesh.template add_property_map<Polygon_mesh::Edge_index, std::set<const Kernel::Plane_3*> >("e:supp_plane").first;//每个边的支撑平面

	typename Polygon_mesh::template Property_map<Polygon_mesh::Vertex_index, std::set<const Kernel::Plane_3*> > vertex_supporting_planes
		= mesh.template add_property_map<Polygon_mesh::Vertex_index, std::set<const Kernel::Plane_3*> >("v:supp_plane").first;//每个顶点的支撑平面

	//为网络中的每个面分配原始平面
	for (auto fd : mesh.faces())
	{
		Polygon_mesh::Halfedge_index h = mesh.halfedge(fd);
		Polygon_mesh::Vertex_index va = mesh.target(h);        const Kernel::Point_3& pa = mesh.points()[va]; h = mesh.next(h);
		Polygon_mesh::Vertex_index vb = mesh.target(h);        const Kernel::Point_3& pb = mesh.points()[vb]; h = mesh.next(h);
		Polygon_mesh::Vertex_index vc = mesh.target(h);        const Kernel::Point_3& pc = mesh.points()[vc];
		const Kernel::Plane_3* plane = new Kernel::Plane_3(pa, pb, pc);//使用三个顶点的坐标构造一个平面
		_supporting_planes.push_back(plane);
		face_supporting_planes[fd] = plane;
	}

	//为网络中的每个边分配原始平面
	for (auto ed : mesh.edges())
	{
		//获取每个边的两个半边
		Polygon_mesh::Halfedge_index h1 = mesh.halfedge(ed);
		Polygon_mesh::Halfedge_index h2 = mesh.opposite(h1);

		//获取半边各自所属的两个平面
		Polygon_mesh::Face_index f1 = mesh.face(h1);
		Polygon_mesh::Face_index f2 = mesh.face(h2);
		CGAL_assertion(f1 != Polygon_mesh::null_face());//bbox网格已关闭
		CGAL_assertion(f2 != Polygon_mesh::null_face());//bbox网格已关闭

		//使用这两个原始平面构造新的平面，存入 edge_supporting_planes 
		const Kernel::Plane_3* plane1 = face_supporting_planes[f1];
		const Kernel::Plane_3* plane2 = face_supporting_planes[f2];
		CGAL_assertion(plane1 && plane2 && plane1 != plane2);

		edge_supporting_planes[ed].insert(plane1);
		edge_supporting_planes[ed].insert(plane2);
		CGAL_assertion(edge_supporting_planes[ed].size() == 2);
	}

	//为网格中的每个顶点分配原始平面
	for (auto vd : mesh.vertices())
	{
		CGAL_assertion(vertex_supporting_planes[vd].size() == 0);
		CGAL::Halfedge_around_target_circulator<Polygon_mesh> hbegin(vd, mesh), done(hbegin);
		do {
			Polygon_mesh::Halfedge_index h = *hbegin;//获取每个顶点的所有半边
			Polygon_mesh::Face_index f = mesh.face(h);//获取该半边所属的面

			//基于该原始平面构造新的平面，存入 vertex_supporting_planes
			const Kernel::Plane_3* plane = face_supporting_planes[f];
			vertex_supporting_planes[vd].insert(plane);
			++hbegin;
		} while (hbegin != done);
		CGAL_assertion(vertex_supporting_planes[vd].size() == 3);
	}

	std::sort(_supporting_planes.begin(), _supporting_planes.end());

	CGAL_assertion(mesh.is_valid());//检查多边形网络，确保其是有效的
}

void HypothesisPlaneProcess::BuildMeshCandidate(Polygon_mesh & bbox_mesh, Polygon_mesh & mesh_candidate)
{
	//1、定义 mesh_bbox 、 mesh_candidate 的属性
	 Polygon_mesh::Property_map<Polygon_mesh::Edge_index, std::set<const Kernel::Plane_3*> > bbox_edge_supporting_planes
		= bbox_mesh.property_map<Polygon_mesh::Edge_index, std::set<const Kernel::Plane_3*> >("e:supp_plane").first;//每个边的所有支持平面

	Polygon_mesh::Property_map<Polygon_mesh::Vertex_index, std::set<const Kernel::Plane_3*> > bbox_vertex_supporting_planes
		= bbox_mesh.property_map<Polygon_mesh::Vertex_index, std::set<const Kernel::Plane_3*> >("v:supp_plane").first;//每个顶点的所有支持平面


	mesh_candidate.clear();//清空 mesh_candidate 网格
	Polygon_mesh::Property_map<Polygon_mesh::Face_index, const Kernel::Plane_3*> face_supporting_planes
		= mesh_candidate.add_property_map<Polygon_mesh::Face_index, const Kernel::Plane_3*>("f:supp_plane").first;//每个面的支持平面

	Polygon_mesh::Property_map<Polygon_mesh::Face_index, CGAL::internal::Planar_segment<Kernel>*> face_supporting_segments
		= mesh_candidate.add_property_map<Polygon_mesh::Face_index, CGAL::internal::Planar_segment<Kernel>*>("f:supp_segment").first;//每个面的支持平面片段

	Polygon_mesh::Property_map<Polygon_mesh::Edge_index, std::set<const Kernel::Plane_3*> > edge_supporting_planes
		= mesh_candidate.add_property_map<Polygon_mesh::Edge_index, std::set<const Kernel::Plane_3*> >("e:supp_plane").first;//每个边的所有支持平面

	Polygon_mesh::Property_map<Polygon_mesh::Vertex_index, std::set<const Kernel::Plane_3*> > vertex_supporting_planes
		= mesh_candidate.add_property_map<Polygon_mesh::Vertex_index, std::set<const Kernel::Plane_3*> >("v:supp_plane").first;//每个顶点的所有支持平面

	//2023.11.23new:添加新的属性映射以表示平面所属的区域
	Polygon_mesh::Property_map<Polygon_mesh::Face_index, std::size_t> face_region_map
		= mesh_candidate.add_property_map<Polygon_mesh::Face_index, std::size_t>("f:region_map").first;//每个面的所属区域


	for (std::size_t i = 0; i < _cloud_plane->planar_segments().size(); ++i)
	{
		CGAL::internal::Planar_segment<Kernel>* g = _cloud_plane->planar_segments()[i];
		const Kernel::Plane_3* cutting_plane = _cloud_plane->planar_segments()[i]->supporting_plane();

		std::vector<Kernel::Point_3> intersecting_points;//容器_交叉点
		std::vector< std::set<const Kernel::Plane_3*> > intersecting_points_source_planes;//容器_交叉点的源平面


		//2、对每个平面，遍历 mesh_bbox 的所有边，判断相交的点并存储
		for (auto ed : bbox_mesh.edges())
		{
			//获取边的两个顶点
			Polygon_mesh::Vertex_index sd = bbox_mesh.vertex(ed, 0);
			Polygon_mesh::Vertex_index td = bbox_mesh.vertex(ed, 1);
			const Kernel::Point_3& s = bbox_mesh.points()[sd];
			const Kernel::Point_3& t = bbox_mesh.points()[td];

			CGAL::Oriented_side ss = cutting_plane->oriented_side(s);//确定内外侧
			CGAL::Oriented_side st = cutting_plane->oriented_side(t);

			//判断平面 plane_cutting 是否与线段 st 相交
			if ((ss == CGAL::ON_POSITIVE_SIDE && st == CGAL::ON_NEGATIVE_SIDE)
				|| (ss == CGAL::ON_NEGATIVE_SIDE && st == CGAL::ON_POSITIVE_SIDE))
			{
				CGAL::Object obj = CGAL::intersection(*cutting_plane, Kernel::Line_3(s, t));//获取相交点
				if (const Kernel::Point_3* p = CGAL::object_cast<Kernel::Point_3>(&obj))
				{
					intersecting_points.push_back(*p);//相交点添加至 intersecting_points
					std::set<const Kernel::Plane_3*> planes = bbox_edge_supporting_planes[ed];//获取与该相交点相关的所有支持平面

					//平面添加至 planes 、 intersecting_points_source_planes
					planes.insert(cutting_plane);
					CGAL_assertion(planes.size() == 3);
					intersecting_points_source_planes.push_back(planes);
				}
				else
					std::cerr << "错误：平面 plane_cutting 与线段 st 不存在交集" << std::endl;
			}
			else
			{
				//plane_cutting 与线段 s-t 不相交
				//判断 plane_cutting 是否与其中一个顶点相交
				if (ss == CGAL::ON_ORIENTED_BOUNDARY && st != CGAL::ON_ORIENTED_BOUNDARY)
				{
					intersecting_points.push_back(s);
					const std::set<const Kernel::Plane_3*>& planes = bbox_vertex_supporting_planes[sd];
					CGAL_assertion(planes.size() == 3);
					intersecting_points_source_planes.push_back(planes);
				}
				else if (st == CGAL::ON_ORIENTED_BOUNDARY && ss != CGAL::ON_ORIENTED_BOUNDARY)
				{
					intersecting_points.push_back(t);
					const std::set<const Kernel::Plane_3*>& planes = bbox_vertex_supporting_planes[td];
					CGAL_assertion(planes.size() == 3);
					intersecting_points_source_planes.push_back(planes);
				}
			}
		}


		//3、确定 intersecting_points 中的凸包点
		if (intersecting_points.size() >= 3)
		{
			//将 intersecting_points 的所有点转为二维坐标，存入 pts
			std::list<Kernel::Point_3> pts;
			for (std::size_t i = 0; i < intersecting_points.size(); ++i)
			{
				const Kernel::Point_3& p = intersecting_points[i];
				const Kernel::Point_2& q = cutting_plane->to_2d(p);//返回仿射变换下的投影的图像点，该投影映射到 XY 平面上，并删除 z 坐标
				pts.push_back(Kernel::Point_3(q.x(), q.y(), FT(i))); // z 分量存储点索引
			}

			//计算 pts 的凸包点集，存入 hull
			std::list<Kernel::Point_3> hull;
			CGAL::convex_hull_2(pts.begin(), pts.end(), std::back_inserter(hull), CGAL::Projection_traits_xy_3<Kernel>());

			//将所有凸包点从 intersecting_points 中存入 point_convexHull 和 ch_source_planes 
			std::vector<Kernel::Point_3> point_convexHull;//容器_凸包点
			std::vector< std::set<const Kernel::Plane_3*> > ch_source_planes;//容器_凸包平面
			for (typename std::list<Kernel::Point_3>::iterator it = hull.begin(); it != hull.end(); ++it)
			{
				std::size_t idx = std::size_t(it->z());
				point_convexHull.push_back(intersecting_points[idx]);
				ch_source_planes.push_back(intersecting_points_source_planes[idx]);
			}


			//4、 mesh_candidate 网络创建新面
			if (point_convexHull.size() >= 3)
			{
				std::vector<Polygon_mesh::Vertex_index> descriptors;
				for (std::size_t j = 0; j < point_convexHull.size(); ++j)
				{
					Polygon_mesh::Vertex_index vd = mesh_candidate.add_vertex(point_convexHull[j]);
					descriptors.push_back(vd);
					vertex_supporting_planes[vd] = ch_source_planes[j];
					CGAL_assertion(vertex_supporting_planes[vd].size() == 3);
				}

				//将该面上的所有顶点的 supporting_planes 属性设置为 ch_source_planes
				Polygon_mesh::Face_index fd = mesh_candidate.add_face(descriptors);
				face_supporting_segments[fd] = _cloud_plane->planar_segments()[i];
				face_supporting_planes[fd] = cutting_plane;

				//new:为该面添加对应的点云区域下标
				face_region_map[fd] = i;
				//std::cout << "i: " << i << "  region_map_i: " << _cloud_plane->property_map<int>("region_map").second[&i] << std::endl;


				//为该面的每个边缘设置 supporting_planes 属性
				CGAL::Halfedge_around_face_circulator<Polygon_mesh> hbegin(mesh_candidate.halfedge(fd), mesh_candidate), done(hbegin);
				do {
					Polygon_mesh::Halfedge_index hd = *hbegin;
					Polygon_mesh::Edge_index ed = mesh_candidate.edge(hd);

					Polygon_mesh::Vertex_index s_vd = mesh_candidate.source(hd);
					Polygon_mesh::Vertex_index t_vd = mesh_candidate.target(hd);
					const std::set<const Kernel::Plane_3*>& s_planes = vertex_supporting_planes[s_vd];
					const std::set<const Kernel::Plane_3*>& t_planes = vertex_supporting_planes[t_vd];
					std::set<const Kernel::Plane_3*> common_planes;
					std::set_intersection(s_planes.begin(), s_planes.end(), t_planes.begin(), t_planes.end(), std::inserter(common_planes, common_planes.begin()));
					if (common_planes.size() == 2)
					{
						//如果两个顶点的 supporting_planes 属性都包含两个平面，则将这两个平面设置为该边缘的 supporting_planes 属性
						CGAL_assertion(edge_supporting_planes[ed].size() == 0);
						edge_supporting_planes[ed] = common_planes;
						CGAL_assertion(edge_supporting_planes[ed].size() == 2);
					}
					else//如果两个顶点的 supporting_planes 属性不都包含两个平面，则说明发生了拓扑错误
						std::cerr << "topological error 拓扑错误" << std::endl;

					++hbegin;
				} while (hbegin != done);
			}
		}

	}//所有的侯选边遍历完毕


}

void HypothesisPlaneProcess::SegMeshCandidate(Polygon_mesh & mesh_candidate)
{
	//1、预先计算平面三元组（依次选取相邻顺序的三个平面）所有潜在的三重交集点（三平面交于一点的点集）
	std::vector<const Kernel::Point_3*> intersecting_points_;//容器_三重交点//8.29_似乎仅做统计，后续用不上
	_triplet_intersections.clear();
	if (_supporting_planes.size() < 4)
		return;//少于4个平面将无法构成闭合网络，退出本函数

	for (std::size_t i = 0; i < _supporting_planes.size(); ++i)
	{
		const Kernel::Plane_3* plane1 = _supporting_planes[i];
		for (std::size_t j = i + 1; j < _supporting_planes.size(); ++j)
		{
			const Kernel::Plane_3* plane2 = _supporting_planes[j];
			for (std::size_t k = j + 1; k < _supporting_planes.size(); ++k)
			{
				const Kernel::Plane_3* plane3 = _supporting_planes[k];
				CGAL_assertion(plane1 < plane2 && plane2 < plane3);//断言_三个平面的指针是否按升序排列的

				if (plane1 == plane2 || plane1 == plane3 || plane2 == plane3)
					continue;//任意两个平面相同，则不会有三重交点，跳过此循环

				//计算三个平面的交集，如果交集的对象可以转换为一个Kernel::Point_3类型的指针，则三个平面相交于一点
				CGAL::Object obj = CGAL::intersection(*plane1, *plane2, *plane3);
				if (const Kernel::Point_3* pt = CGAL::object_cast<Kernel::Point_3>(&obj))
				{
					//复制三重交点的坐标，存储点
					Kernel::Point_3* new_point = new Kernel::Point_3(*pt);
					_triplet_intersections[plane1][plane2][plane3] = new_point;
					intersecting_points_.push_back(new_point);
				}
				else
					continue;//未得到交点的两种原因：（1）面是平行的；（2）面相交于同一条线上。

			}
		}
	}
	///输出测试_三重交点
	//std::cout << "交点个数：" << intersecting_points_.size() << std::endl;
	//for (std::vector<const Kernel::Point_3*>::iterator it = intersecting_points_.begin(); it != intersecting_points_.end(); it++)
	//	std::cout << (*it)->x() << "  "
	//		<< (*it)->y() << "  "
	//		<< (*it)->z() << "  " << std::endl;


	//2、求取候选平面之间的相交关系（确定相交面的各个交点）
	//不能直接使用 const Polygon_mesh::Face_range& all_faces = mesh.faces();
	std::vector<Polygon_mesh::Face_index> all_faces(mesh_candidate.faces().begin(), mesh_candidate.faces().end());//容器_所有候选面

	typename Polygon_mesh::template Property_map<Polygon_mesh::Face_index, const Kernel::Plane_3*> face_supporting_planes
		= mesh_candidate.template property_map<Polygon_mesh::Face_index, const Kernel::Plane_3*>("f:supp_plane").first;//属性映射_每个平面的支撑面

	for (std::size_t i = 0; i < all_faces.size(); ++i)
	{
		//获取第i个面的索引和所在平面
		Polygon_mesh::Face_index index_face_i = all_faces[i];
		const Kernel::Plane_3* plane_face_i = face_supporting_planes[index_face_i];

		//获取相应的属性映射
		typename Polygon_mesh::template Property_map<Polygon_mesh::Face_index, const Kernel::Plane_3*> face_supporting_planes =
			mesh_candidate.template property_map<Polygon_mesh::Face_index, const Kernel::Plane_3*>("f:supp_plane").first;//属性映射_每个面的支撑平面

		typename Polygon_mesh::template Property_map<Polygon_mesh::Face_index, CGAL::internal::Planar_segment<Kernel>*> face_supporting_segments =
			mesh_candidate.template property_map<Polygon_mesh::Face_index, CGAL::internal::Planar_segment<Kernel>*>("f:supp_segment").first;//属性映射_每个面的支撑平面段


		//2.1、遍历所有候选面，将与第i个面相交的面存入集合 intersecting_faces
		std::set<Polygon_mesh::Face_index> intersecting_faces;
		for (auto f : mesh_candidate.faces())
		{
			//如果f与第i个面相同/f与第i个面在同一平面段上/f与第i个面在同一平面上，则两者不会相交，跳过
			if (f == index_face_i ||
				face_supporting_segments[f] == face_supporting_segments[index_face_i] ||
				face_supporting_planes[f] == face_supporting_planes[index_face_i])
				continue;

			//2.2、判断面 f 是否与平面 plane_cutting 相交
			const Kernel::Plane_3* plane_cutting = face_supporting_planes[index_face_i];//获取给定面作为切割面
			CGAL_assertion(plane_cutting != nullptr);//检查平面是否为空指针

			std::vector<Polygon_mesh::Vertex_index> existing_vts;//平面相交的已有顶点下标
			std::vector<EdgePos> new_vts;//平面相交的新生成的顶点
			compute_intersections(mesh_candidate, f, plane_cutting, existing_vts, new_vts);//计算交点

			bool is_intersect_fPlane = false;
			if (existing_vts.size() == 2)
			{
				if (!halfedge_exists(existing_vts[0], existing_vts[1], mesh_candidate))
					is_intersect_fPlane = true;//2个已有的顶点与plane_cutting相交且两点不在同一条边上
			}
			else if (existing_vts.size() + new_vts.size() == 2)
				is_intersect_fPlane = true;//已有+新生成的共两个点与plane_cutting相交

			if (is_intersect_fPlane)
				intersecting_faces.insert(f);//如果f和第i个面所在的平面相交,将f加入到intersecting_faces

			///输出测试_面的交点-------
			//if (existing_vts.size() != 0)
			//{
			//	for (int i = 0; i < existing_vts.size(); ++i)
			//	{
			//		Kernel::Point_3 temp = mesh_candidate.points()[existing_vts[i]];
			//		std::cout << "existing_vts" << temp.x() << "  " << temp.y() << "  " << temp.z() << std::endl;
			//		//point_intersection_test.push_back(temp);
			//	}
			//}
			//if (new_vts.size() != 0)
			//{
			//	for (int i = 0; i < new_vts.size(); ++i)
			//	{
			//		Kernel::Point_3 temp(new_vts[i].pos->x(), new_vts[i].pos->y(), new_vts[i].pos->z());
			//		std::cout << "existing_vts" << temp.x() << "  " << temp.y() << "  " << temp.z() << std::endl;
			//		//point_intersection_test.push_back(temp);
			//	}
			//}
			//std::cout << std::endl;
			///------------------------

		}

		if (intersecting_faces.empty())
			continue;//对第i个候选面，没有与之相交的候选面，则跳过
		std::vector<Polygon_mesh::Face_index> cutting_faces(intersecting_faces.begin(), intersecting_faces.end());//容器_与第i个候选面相交的其他候选面

		///输出测试_相交的面
		//std::cout << "相交平面：" << cutting_faces.size() << std::endl;
		//for (std::vector<Polygon_mesh::Face_index>::iterator it = cutting_faces.begin(); it != cutting_faces.end(); ++it)
		//	std::cout << (*it) << std::endl;


		//3、使用与 index_face_i 相交的平面分割 index_face_i 
		//每次切割后，原来的面不再存在，它被多个块取代，每个块都将被另一个平面切割
		std::vector<Polygon_mesh::Face_index> faces_to_be_cut;//存储需要被分割的面
		faces_to_be_cut.push_back(index_face_i);//先添加当前的第i个面
		while (!intersecting_faces.empty())
		{
			//以第一个相交面作为切割面
			Polygon_mesh::Face_index cutting_face = *(intersecting_faces.begin());//切割面的索引
			const Kernel::Plane_3* cutting_plane = face_supporting_planes[cutting_face];//获取切割面所在平面

			std::set<Polygon_mesh::Face_index> new_faces;                //存储切割后产生的新面
			std::set<Polygon_mesh::Face_index> remained_faces;        //存储还需被切割的面
			for (std::size_t j = 0; j < faces_to_be_cut.size(); ++j)
			{
				Polygon_mesh::Face_index current_face = faces_to_be_cut[j];//被切割的面
				std::vector<Polygon_mesh::Face_index> tmp = split_plane(current_face, cutting_plane, mesh_candidate);//调用split_plane函数对这个面进行切割
				new_faces.insert(tmp.begin(), tmp.end());//把这些新的面加入到new_faces中
				if (tmp.empty())
				{
					remained_faces.insert(current_face);//没有发生切割，将该面加入到 remained_faces 待后续进行切割
				}
			}

			//把 new_faces 中的面作为下一轮需要被切割的面
			faces_to_be_cut = std::vector<Polygon_mesh::Face_index>(new_faces.begin(), new_faces.end());

			//把 remained_faces 中的面也加入到需要被切割的面
			faces_to_be_cut.insert(faces_to_be_cut.end(), remained_faces.begin(), remained_faces.end());

			intersecting_faces.erase(cutting_face);//将被切割面从相交的面中移除
		}


		//4、所有和第i个面相交的面 cutting_faces 都会被第i个面所在的平面 index_face_i 切割
		for (std::size_t j = 0; j < cutting_faces.size(); ++j)
			split_plane(cutting_faces[j], plane_face_i, mesh_candidate);

	}

	CGAL_assertion(mesh_candidate.is_valid());
}

void HypothesisPlaneProcess::SimpleMeshCandidate(Polygon_mesh & mesh_candidate)
{
	std::vector<Polygon_mesh::Face_index> faces_covered;//点云覆盖区域及邻近的面下标
	for (size_t i_segment = 0; i_segment < _cloud_plane->planar_segments().size(); i_segment++)
	{
		//1、获取区域点的覆盖区域
		CGAL::internal::Planar_segment<Kernel>* s = _cloud_plane->planar_segments()[i_segment];//指针_分段区域

		//获取区域点云
		std::vector<Kernel::Point_3> points_region;
		for (std::size_t i = 0; i < s->size(); ++i)
		{
			std::size_t idx = s->at(i);//点索引
			Kernel::Point_3 p = _cloud_plane->point_map()[idx];
			points_region.push_back(p);
		}

		//三维点投影到平面
		const Kernel::Plane_3* plane_fitting = _cloud_plane->planar_segments()[i_segment]->supporting_plane();//拟合平面
		std::vector<Kernel::Point_2> points_region_2;
		for (const Kernel::Point_3& p : points_region)
			points_region_2.push_back(plane_fitting->to_2d(p));


		//2、遍历多边形上的每个面，将包含投影点及相关的面标记到 faces_covered
		std::vector<Polygon_mesh::Face_index> faces_polygon;//存储多边形平面下标
		for (const Polygon_mesh::Face_index& f : mesh_candidate.faces())
			if (mesh_candidate.property_map<Polygon_mesh::Face_index, std::size_t>("f:region_map").first[f] == i_segment)
				faces_polygon.push_back(f);


		GeometryCal obj_gc;
		std::vector< Polygon_mesh::Vertex_index> vi_covered;
		for (const Polygon_mesh::Face_index& f : faces_polygon)
		{
			std::vector<Kernel::Point_3> points_polygon;
			std::vector<Polygon_mesh::Vertex_index> vi_polygon;

			//类模板，创建顶点环绕面的两个索引器
			//vcirc (h, sm)表示从h开始环绕mesh中面的顶点；done(vcirc)用于判断循环是否结束
			CGAL::Vertex_around_face_circulator<Polygon_mesh> vcirc(mesh_candidate.halfedge(f), mesh_candidate), done(vcirc);
			do {
				Polygon_mesh::Vertex_index vi = *vcirc;//顶点索引
				Kernel::Point_3 p = mesh_candidate.point(vi);//顶点坐标
				vi_polygon.push_back(vi);
				points_polygon.push_back(p);

				++vcirc;
			} while (vcirc != done);

			//多边形顶点投影到平面
			std::vector<Kernel::Point_2> points_polygon_2;
			for (const Kernel::Point_3& p : points_polygon)
				points_polygon_2.push_back(plane_fitting->to_2d(p));

			//判断是否有区域点云落在多边形内
			for (Kernel::Point_2 p : points_region_2)
				if (obj_gc.isPointInPolygon(p, points_polygon_2))
				{
					//遍历该面 f 的顶点，存入
					for (const Polygon_mesh::Vertex_index& vi : vi_polygon)
						if (std::find(vi_covered.begin(), vi_covered.end(), vi) == vi_covered.end())
							vi_covered.push_back(vi);

					break;//该面 f 的顶点数据已存储，无需再遍历 points_region_2
				}
		}


		for (int i = 0; i < vi_covered.size(); ++i)
			for (const Polygon_mesh::Face_index& f : faces_polygon)
			{
				if (std::find(faces_covered.begin(), faces_covered.end(), f) != faces_covered.end())
					continue;//该面 f 已存入 faces_covered，跳过

				CGAL::Vertex_around_face_circulator<Polygon_mesh> vcirc(mesh_candidate.halfedge(f), mesh_candidate), done(vcirc);
				do {
					Polygon_mesh::Vertex_index vi = *vcirc;//顶点索引
					if (vi == vi_covered[i])
					{
						faces_covered.push_back(f);//faces_covered 中未存在面 f ，将其存入
						break;
					}

					++vcirc;
				} while (vcirc != done);
			}
	}


	//3、统计每个顶点的使用次数，剔除边缘面
	bool isRemove = false;//是否有平面段被合并
	do {
		isRemove = false;

		//统计所有f中v出现的次数
		std::unordered_map<Kernel::Point_3, int> map_vNums;
		for (const Polygon_mesh::Face_index& f : faces_covered)
		{
			CGAL::Vertex_around_face_circulator<Polygon_mesh> vcirc(mesh_candidate.halfedge(f), mesh_candidate), done(vcirc);
			do {
				Polygon_mesh::Vertex_index vi = *vcirc;//顶点索引
				Kernel::Point_3 p = mesh_candidate.point(vi);//顶点坐标

				map_vNums[p]++;

				++vcirc;
			} while (vcirc != done);
		}

		///输出测试
		//for (const auto& pair : map_vNums)
		//	std::cout << pair.first << " 出现了 " << pair.second << " 次" << std::endl;

		//登记有边缘面
		std::vector<Kernel::Point_3> vi_remove;
		for (auto pair : map_vNums)
			if (pair.second == 1)
				vi_remove.push_back(pair.first);

		///输出测试
		//std::cout << "vi_remove:" << vi_remove.size() << std::endl
		//	<< "map_vNums:" << map_vNums.size() << std::endl << std::endl;

		//存在需要被剔除的面
		if (!vi_remove.empty())
		{
			//记录需要被剔除的f
			std::vector<Polygon_mesh::Face_index> faces_remove;
			for (int i = 0; i < vi_remove.size(); ++i)
				for (const Polygon_mesh::Face_index& f : faces_covered)
				{
					if (std::find(faces_remove.begin(), faces_remove.end(), f) != faces_remove.end())
						continue;//该面 f 已存入 faces_remove，跳过

					CGAL::Vertex_around_face_circulator<Polygon_mesh> vcirc(mesh_candidate.halfedge(f), mesh_candidate), done(vcirc);
					do {
						Polygon_mesh::Vertex_index vi = *vcirc;//顶点索引
						Kernel::Point_3 p = mesh_candidate.point(vi);//顶点坐标

						if (p == vi_remove[i])
						{
							faces_remove.push_back(f);//faces_remove 中未存在面 f ，将其存入
							break;
						}

						++vcirc;
					} while (vcirc != done);
				}

			//移除faces_covered中对应的f
			for (const Polygon_mesh::Face_index& f : faces_remove)
				faces_covered.erase(std::remove(faces_covered.begin(), faces_covered.end(), f), faces_covered.end());

			isRemove = true;
		}

	} while (isRemove);//isRemove 若为false，则此时所有平面都不适合被移除，跳出循环


	//4、简化面
	std::vector<Polygon_mesh::Face_index> faces_delete;//面索引_记录需要删除的面
	for (const Polygon_mesh::Face_index& f : mesh_candidate.faces())
	{
		if (std::find(faces_covered.begin(), faces_covered.end(), f) != faces_covered.end())
			continue;//该面f属于faces_covered，不用删除

		faces_delete.push_back(f);
	}

	for (std::size_t i = 0; i < faces_delete.size(); ++i)
	{
		Polygon_mesh::Face_index f = faces_delete[i];
		Polygon_mesh::Halfedge_index h = mesh_candidate.halfedge(f);
		CGAL::Euler::remove_face(h, mesh_candidate);//根据半边索引删除相应面
	}

	mesh_candidate.collect_garbage();//重新组织网络，移除被删除的面

	///输出测试
	//for (const Polygon_mesh::Face_index& f : mesh_inside.faces())
	//	std::cout << *mesh_inside.property_map<Polygon_mesh::Face_index, const Kernel::Plane_3*>("f:supp_plane").first[f] << "  "
	//	<< f << std::endl;


#if 0

	//似乎无效
	_triplet_intersections.clear();
	std::vector<const Kernel::Plane_3*> plane_support_new;
	for (const Polygon_mesh::Face_index& f : mesh_candidate.faces())
	{
		const Kernel::Plane_3* supp_plane =
			mesh_candidate.property_map<Polygon_mesh::Face_index, const Kernel::Plane_3*>("f:supp_plane").first[f];

		if (std::find(plane_support_new.begin(), plane_support_new.end(), supp_plane) == plane_support_new.end())
			plane_support_new.push_back(supp_plane);
	}
	std::sort(plane_support_new.begin(), plane_support_new.end());

	//预先计算平面三元组（依次选取相邻顺序的三个平面）所有潜在的三重交集点（三平面交于一点的点集）
	std::vector<const Kernel::Point_3*> intersecting_points_;//容器_三重交点//8.29_似乎仅做统计，后续用不上
	for (std::size_t i = 0; i < plane_support_new.size(); ++i)
	{
		const Kernel::Plane_3* plane1 = plane_support_new[i];
		for (std::size_t j = i + 1; j < plane_support_new.size(); ++j)
		{
			const Kernel::Plane_3* plane2 = plane_support_new[j];
			for (std::size_t k = j + 1; k < plane_support_new.size(); ++k)
			{
				const Kernel::Plane_3* plane3 = plane_support_new[k];
				CGAL_assertion(plane1 < plane2 && plane2 < plane3);//断言_三个平面的指针是否按升序排列的

				if (plane1 == plane2 || plane1 == plane3 || plane2 == plane3)
					continue;//任意两个平面相同，则不会有三重交点，跳过此循环

				//计算三个平面的交集，如果交集的对象可以转换为一个Kernel::Point_3类型的指针，则三个平面相交于一点
				CGAL::Object obj = CGAL::intersection(*plane1, *plane2, *plane3);
				if (const Kernel::Point_3* pt = CGAL::object_cast<Kernel::Point_3>(&obj))
				{
					//复制三重交点的坐标，存储点
					Kernel::Point_3* new_point = new Kernel::Point_3(*pt);
					_triplet_intersections[plane1][plane2][plane3] = new_point;
					intersecting_points_.push_back(new_point);
				}
				else
					continue;//未得到交点的两种原因：（1）面是平行的；（2）面相交于同一条线上。

			}
		}
	}
	///输出测试_三重交点
	std::cout << "交点个数：" << intersecting_points_.size() << std::endl;

#endif // 重新计算潜在的三重交点

}

void HypothesisPlaneProcess::SimpleMeshCandidate_V2(Polygon_mesh & mesh_candidate)
{
	std::vector< Kernel::Point_3> vi_covered;//点云覆盖区域及邻近的点下标
	for (size_t i_segment = 0; i_segment < _cloud_plane->planar_segments().size(); i_segment++)
	{
		//1、获取区域点的覆盖区域
		CGAL::internal::Planar_segment<Kernel>* s = _cloud_plane->planar_segments()[i_segment];//指针_分段区域

		//获取区域点云
		std::vector<Kernel::Point_3> points_region;
		for (std::size_t i = 0; i < s->size(); ++i)
		{
			std::size_t idx = s->at(i);//点索引
			Kernel::Point_3 p = _cloud_plane->point_map()[idx];
			points_region.push_back(p);
		}

		//三维点投影到平面
		const Kernel::Plane_3* plane_fitting = _cloud_plane->planar_segments()[i_segment]->supporting_plane();//拟合平面
		std::vector<Kernel::Point_2> points_region_2;
		for (const Kernel::Point_3& p : points_region)
			points_region_2.push_back(plane_fitting->to_2d(p));


		//2、遍历多边形上的每个面，将包含投影点及相关的面标记到 faces_covered
		std::vector<Polygon_mesh::Face_index> faces_polygon;//存储多边形平面下标
		for (const Polygon_mesh::Face_index& f : mesh_candidate.faces())
			if (mesh_candidate.property_map<Polygon_mesh::Face_index, std::size_t>("f:region_map").first[f] == i_segment)
				faces_polygon.push_back(f);


		GeometryCal obj_gc;
		for (const Polygon_mesh::Face_index& f : faces_polygon)
		{
			std::vector<Kernel::Point_3> points_polygon;
			std::vector<Polygon_mesh::Vertex_index> vi_polygon;

			//类模板，创建顶点环绕面的两个索引器
			//vcirc (h, sm)表示从h开始环绕mesh中面的顶点；done(vcirc)用于判断循环是否结束
			CGAL::Vertex_around_face_circulator<Polygon_mesh> vcirc(mesh_candidate.halfedge(f), mesh_candidate), done(vcirc);
			do {
				Polygon_mesh::Vertex_index vi = *vcirc;//顶点索引
				Kernel::Point_3 p = mesh_candidate.point(vi);//顶点坐标
				vi_polygon.push_back(vi);
				points_polygon.push_back(p);

				++vcirc;
			} while (vcirc != done);

			//多边形顶点投影到平面
			std::vector<Kernel::Point_2> points_polygon_2;
			for (const Kernel::Point_3& p : points_polygon)
				points_polygon_2.push_back(plane_fitting->to_2d(p));

			//判断是否有区域点云落在多边形内
			for (Kernel::Point_2 p : points_region_2)
				if (obj_gc.isPointInPolygon(p, points_polygon_2))
				{
					//遍历该面 f 的顶点，存入
					for (const Polygon_mesh::Vertex_index& vi : vi_polygon)
					{
						Kernel::Point_3 p = mesh_candidate.point(vi);//顶点坐标
						if (std::find(vi_covered.begin(), vi_covered.end(), p) == vi_covered.end())
							vi_covered.push_back(p);
					}

					break;//该面 f 的顶点数据已存储，无需再遍历 points_region_2
				}
		}
	}


	//3、统计覆盖区域的顶点并保留关联的面，剔除边缘面
	std::vector<Polygon_mesh::Face_index> faces_covered;//点云覆盖区域及邻近的面下标
	for (const Polygon_mesh::Face_index& f : mesh_candidate.faces())
	{
		if (std::find(faces_covered.begin(), faces_covered.end(), f) != faces_covered.end())
			continue;//该面 f 已存入 faces_covered，跳过

		CGAL::Vertex_around_face_circulator<Polygon_mesh> vcirc(mesh_candidate.halfedge(f), mesh_candidate), done(vcirc);
		do {
			Polygon_mesh::Vertex_index vi = *vcirc;//顶点索引
			Kernel::Point_3 p = mesh_candidate.point(vi);//顶点坐标

			if (std::find(vi_covered.begin(), vi_covered.end(), p) != vi_covered.end())
			{
				faces_covered.push_back(f);//faces_covered 中未存在面 f ，将其存入
				break;
			}

			++vcirc;
		} while (vcirc != done);
	}

	///输出测试
	//faces_covered.clear();
	//for (const Polygon_mesh::Face_index& f : mesh_candidate.faces())
	//{
	//	if (f.idx() == 235 || f.idx() == 444)
	//		continue;
	//	faces_covered.push_back(f);
	//}


	//4、简化面
	std::vector<Polygon_mesh::Face_index> faces_delete;//面索引_记录需要删除的面
	for (const Polygon_mesh::Face_index& f : mesh_candidate.faces())
	{
		if (std::find(faces_covered.begin(), faces_covered.end(), f) != faces_covered.end())
			continue;//该面f属于faces_covered，不用删除

		faces_delete.push_back(f);
	}

	for (std::size_t i = 0; i < faces_delete.size(); ++i)
	{
		Polygon_mesh::Face_index f = faces_delete[i];
		Polygon_mesh::Halfedge_index h = mesh_candidate.halfedge(f);
		CGAL::Euler::remove_face(h, mesh_candidate);//根据半边索引删除相应面
	}

	mesh_candidate.collect_garbage();//重新组织网络，移除被删除的面

	///输出测试
	//for (const Polygon_mesh::Face_index& f : mesh_inside.faces())
	//	std::cout << *mesh_inside.property_map<Polygon_mesh::Face_index, const Kernel::Plane_3*>("f:supp_plane").first[f] << "  "
	//	<< f << std::endl;

#if 0

	//似乎无效
	_triplet_intersections.clear();
	std::vector<const Kernel::Plane_3*> plane_support_new;
	for (const Polygon_mesh::Face_index& f : mesh_candidate.faces())
	{
		const Kernel::Plane_3* supp_plane =
			mesh_candidate.property_map<Polygon_mesh::Face_index, const Kernel::Plane_3*>("f:supp_plane").first[f];

		if (std::find(plane_support_new.begin(), plane_support_new.end(), supp_plane) == plane_support_new.end())
			plane_support_new.push_back(supp_plane);
	}
	std::sort(plane_support_new.begin(), plane_support_new.end());

	//预先计算平面三元组（依次选取相邻顺序的三个平面）所有潜在的三重交集点（三平面交于一点的点集）
	std::vector<const Kernel::Point_3*> intersecting_points_;//容器_三重交点//8.29_似乎仅做统计，后续用不上
	for (std::size_t i = 0; i < plane_support_new.size(); ++i)
	{
		const Kernel::Plane_3* plane1 = plane_support_new[i];
		for (std::size_t j = i + 1; j < plane_support_new.size(); ++j)
		{
			const Kernel::Plane_3* plane2 = plane_support_new[j];
			for (std::size_t k = j + 1; k < plane_support_new.size(); ++k)
			{
				const Kernel::Plane_3* plane3 = plane_support_new[k];
				CGAL_assertion(plane1 < plane2 && plane2 < plane3);//断言_三个平面的指针是否按升序排列的

				if (plane1 == plane2 || plane1 == plane3 || plane2 == plane3)
					continue;//任意两个平面相同，则不会有三重交点，跳过此循环

				//计算三个平面的交集，如果交集的对象可以转换为一个Kernel::Point_3类型的指针，则三个平面相交于一点
				CGAL::Object obj = CGAL::intersection(*plane1, *plane2, *plane3);
				if (const Kernel::Point_3* pt = CGAL::object_cast<Kernel::Point_3>(&obj))
				{
					//复制三重交点的坐标，存储点
					Kernel::Point_3* new_point = new Kernel::Point_3(*pt);
					_triplet_intersections[plane1][plane2][plane3] = new_point;
					intersecting_points_.push_back(new_point);
				}
				else
					continue;//未得到交点的两种原因：（1）面是平行的；（2）面相交于同一条线上。

			}
		}
	}
	///输出测试_三重交点
	std::cout << "交点个数：" << intersecting_points_.size() << std::endl;

#endif // 重新计算潜在的三重交点

}

void HypothesisPlaneProcess::AddMeshProperty(Polygon_mesh & mesh_candidate)
{
	typename Polygon_mesh::template Property_map<Polygon_mesh::Face_index, double> face_roof =
		mesh_candidate.template add_property_map<Polygon_mesh::Face_index, double>("f:dist_boxTop").first;//属性映射_面片属顶面

	//候选面顶点集
	std::vector<Kernel::Point_3> vertices(mesh_candidate.number_of_vertices());
	for (int i = 0; i < mesh_candidate.vertices().size(); ++i)
	{
		const Polygon_mesh::Vertex_index v = mesh_candidate.vertices().first[i];
		vertices[i] = mesh_candidate.points()[v];
	}

	//计算包围盒
	const Kernel::Iso_cuboid_3& box = CGAL::bounding_box(vertices.begin(), vertices.end());//计算 vertices 点集的最小包围盒
	FT dx = box.xmax() - box.xmin();//包围盒xyz方向上的长度
	FT dy = box.ymax() - box.ymin();
	FT dz = box.zmax() - box.zmin();

	//计算包围盒顶面
	std::vector<Kernel::Point_3> point_boxTop;//包围盒顶面点
	point_boxTop.push_back(Kernel::Point_3(box.xmin(), box.ymin(), box.zmax()));
	point_boxTop.push_back(Kernel::Point_3(box.xmax(), box.ymin(), box.zmax()));
	point_boxTop.push_back(Kernel::Point_3(box.xmax(), box.ymax(), box.zmax()));
	point_boxTop.push_back(Kernel::Point_3(box.xmin(), box.ymax(), box.zmax()));
	Kernel::Plane_3 plane_boxTop;
	CGAL::linear_least_squares_fitting_3(point_boxTop.begin(), point_boxTop.end(), plane_boxTop, CGAL::Dimension_tag<0>());

	//遍历候选面集的每个面片
	for (const Polygon_mesh::Face_index& f : mesh_candidate.faces())
	{
		double x_center = 0.0f, y_center = 0.0f, z_center = 0.0f;
		int num_p = 0;

		//累计该面上的点坐标
		CGAL::Vertex_around_face_circulator<Polygon_mesh> vcirc(mesh_candidate.halfedge(f), mesh_candidate), done(vcirc);
		do {
			Polygon_mesh::Vertex_index vi = *vcirc;//顶点索引
			Kernel::Point_3 p = mesh_candidate.point(vi);//顶点坐标

			x_center += p.x();
			y_center += p.y();
			z_center += p.z();
			++num_p;

			++vcirc;
		} while (vcirc != done);

		//面片中心点
		Kernel::Point_3 p_center(x_center / num_p, y_center / num_p, z_center / num_p);

		//点到平面的距离 / 包围盒高度
		//除以包围盒高度是想归一化处理，让值处于[0,1]
		double dist_p2plane = std::sqrt(CGAL::squared_distance(p_center, plane_boxTop)) /(box.zmax() - box.zmin());

		//属性映射赋值
		face_roof[f] = dist_p2plane;
	}

}

void HypothesisPlaneProcess::BuildHypothesisPlane(const PointCloud & cloud_region, CGAL::Surface_mesh<Kernel::Point_3>& mesh_candidate)
{
	CGAL::Timer timer;//时间计数的实例
	timer.start();

	//0、转格式，cloud_region -> cloud_plane
	std::vector<PNI> point_pni;
	for (int i = 0; i < cloud_region.size(); ++i)
	{
		PNI temp;
		temp.get<0>() = cloud_region.point(i);
		temp.get<1>() = cloud_region.normal(i);
		temp.get<2>() = cloud_region.property_map<int>("region_map").first[i];
		point_pni.push_back(temp);
	}

	//2023.8.21 注意，cloud_plane要全程引用（等同于全局变量），不然容易触发~Point_set_with_planes()从而清掉 planar_segments
	CGAL::internal::Point_set_with_planes<Kernel> cloud_plane(point_pni,//点云对象
		CGAL::Nth_of_tuple_property_map<0, PNI>(),//指定类型为Point_map
		CGAL::Nth_of_tuple_property_map<1, PNI>(),//Normal_map
		CGAL::Nth_of_tuple_property_map<2, PNI>());//Plane_index_map

	if (cloud_plane.planar_segments().size() < 4)
	{
		std::cout << "错误：重建闭合表面模型至少需要4个平面聚类，该点云只有 " + std::to_string(cloud_plane.planar_segments().size()) + " 个聚类" << std::endl;
		return;//检查点集的平面段数量，如果小于4，则返回错误消息
	}
	_cloud_plane = &cloud_plane;


	//1、优化点云的平面聚类（通过合并压低数量）
	//std::cout <<"当前点云的平面段数量："<< _cloud_plane->planar_segments().size() << std::endl;
	PlaneRefine();
	//std::cout << "优化后的点云的平面段数量：" << _cloud_plane->planar_segments().size() << std::endl;


	//2、构建点云的包围盒多边形
	CGAL::Surface_mesh<Kernel::Point_3> mesh_bbox;
	BuildMeshBbox(mesh_bbox);


	//3、构建点云的候选表面多边形（跨幅延伸至包围盒的平面）
	BuildMeshCandidate(mesh_bbox, mesh_candidate);


	//4、候选平面相交细分
	//2023.11.23 在包围盒平面与细分的平面上，对面的属性增加分区属性
	SegMeshCandidate(mesh_candidate);
	

	//5.1、候选面集简化处理
	//SimpleMeshCandidate(mesh_candidate);//激进简化，有时重建效果巨差
	SimpleMeshCandidate_V2(mesh_candidate);//保守简化

	//5.2、附加顶面偏好属性
	AddMeshProperty(mesh_candidate);


	//6、计算候选面的置信度
	//该操作为 mesh_candidate 添加了3个属性映射,用于后续的MIP处理：
	//平面的支撑点数	f:num_supporting_points
	//平面面积			f:face_area
	//覆盖面积			f:covered_area
	CGAL::internal::Candidate_confidences<Kernel> conf;
	conf.compute(cloud_plane, mesh_candidate);

	///输出测试_显示新加的三种属性映射
	//std::cout << mesh_candidate.faces().size() << std::endl;
	//for (auto f : mesh_candidate.faces())
	//{
	//	std::cout << mesh_candidate.property_map<Polygon_mesh::Face_index, std::size_t>("f:num_supporting_points").first[f] << std::endl;
	//	std::cout << mesh_candidate.property_map<Polygon_mesh::Face_index, FT>("f:face_area").first[f] << std::endl;
	//	std::cout << mesh_candidate.property_map<Polygon_mesh::Face_index, FT>("f:covered_area").first[f] << std::endl;
	//	std::cout << std::endl;
	//}


	timer.stop();
	std::cout << "有 " << mesh_candidate.faces().size() << " 个候选表面已生成，运行耗时 " << timer.time() << " 秒" << std::endl;


	///输出测试_候选面的颜色属性//弃用
	//for (auto f : mesh_candidate.faces())
	//	std::cout << mesh_candidate.property_map<Polygon_mesh::Face_index, CGAL::IO::Color>("f:color").first[f] << std::endl;
}

void HypothesisPlaneProcess::Coordination2Original(Polygon_mesh & model, double T_x, double T_y, double T_z)
{
	//复用
	GeometryCal obj_geo;
	obj_geo.Coordination2Original(model, T_x, T_y, T_z);//网络平移至原始坐标
}







std::size_t HypothesisPlaneProcess::calNumber_points_on_plane(const CGAL::internal::Planar_segment<Kernel>* s, const Kernel::Plane_3 * plane, FT dist_threshold)
{
	CGAL_assertion(const_cast<CGAL::internal::Planar_segment<Kernel>*>(s)->point_set() == _cloud_plane);//检查参数有效性

	//1、初始化变量
	//获取 point_set_ 对象的 point_map() 成员变量的引用
	//可访问并修改 point_set_ 对象的所有点
	const CGAL::internal::Point_set_with_planes<Kernel>::Point_map& points = _cloud_plane->point_map();

	//2、遍历平面段 s 上的所有点
	std::size_t count = 0;
	for (std::size_t i = 0; i < s->size(); ++i)
	{
		std::size_t idx = s->at(i);
		const Kernel::Point_3& p = points[idx];

		//3、计算点 p 到平面 plane_cutting 的距离
		FT sdist = CGAL::squared_distance(*plane, p);
		FT dist = std::sqrt(sdist);
		if (dist < dist_threshold)//如果距离小于阈值，计数+1
			++count;
	}
	return count;
}

void HypothesisPlaneProcess::plane_merge(CGAL::internal::Planar_segment<Kernel>* s1, CGAL::internal::Planar_segment<Kernel>* s2)
{
	//检查参数有效性
	CGAL_assertion(const_cast<CGAL::internal::Planar_segment<Kernel>*>(s1)->point_set() == _cloud_plane);
	CGAL_assertion(const_cast<CGAL::internal::Planar_segment<Kernel>*>(s2)->point_set() == _cloud_plane);
	//获取 point_set_ 对象的 planar_segments() 成员变量的引用
	//可访问并修改 point_set_ 对象的所有线段
	std::vector< CGAL::internal::Planar_segment<Kernel>* >& segments = _cloud_plane->planar_segments();

	//1、获取平面 s1 和 s2 的所有点的索引
	std::vector<std::size_t> points_indices;
	points_indices.insert(points_indices.end(), s1->begin(), s1->end());
	points_indices.insert(points_indices.end(), s2->begin(), s2->end());

	//2、创建一个新的平面 s，并设置它的所有点为 s1 + s2
	CGAL::internal::Planar_segment<Kernel>* s = new CGAL::internal::Planar_segment<Kernel>(_cloud_plane);
	s->insert(s->end(), points_indices.begin(), points_indices.end());
	s->fit_supporting_plane();//计算线段的支撑平面
	segments.push_back(s);//将新平面 s 添加到平面段 segments

	//3、从点集内删除平面 s1 和 s2
	typename std::vector< CGAL::internal::Planar_segment<Kernel>* >::iterator pos = std::find(segments.begin(), segments.end(), s1);
	if (pos != segments.end())
	{
		CGAL::internal::Planar_segment<Kernel>* tmp = *pos;
		const Kernel::Plane_3* plane = tmp->supporting_plane();
		segments.erase(pos);
		delete tmp;
		delete plane;
	}
	else
		std::cerr << "Fatal error: should not reach here" << std::endl;

	pos = std::find(segments.begin(), segments.end(), s2);
	if (pos != segments.end())
	{
		CGAL::internal::Planar_segment<Kernel>* tmp = *pos;
		const Kernel::Plane_3* plane = tmp->supporting_plane();
		segments.erase(pos);
		delete tmp;
		delete plane;
	}
	else
		std::cerr << "Fatal error: should not reach here" << std::endl;

}

void HypothesisPlaneProcess::compute_intersections(const Polygon_mesh & mesh, Polygon_mesh::Face_index face, const Kernel::Plane_3 * plane_cutting, std::vector<Polygon_mesh::Vertex_index>& existing_vts, std::vector<EdgePos>& new_vts)
{
	existing_vts.clear();
	new_vts.clear();

	typename Polygon_mesh::template Property_map<Polygon_mesh::Face_index, const Kernel::Plane_3*> face_supporting_planes =
		mesh.template property_map<Polygon_mesh::Face_index, const Kernel::Plane_3*>("f:supp_plane").first;//属性映射_每个给定面所在的支撑面
	const Kernel::Plane_3* supporting_plane = face_supporting_planes[face];
	if (supporting_plane == plane_cutting)
		return;//给定面与切割平面相同，则不需计算交点，退出函数

	typename Polygon_mesh::template Property_map<Polygon_mesh::Edge_index, std::set<const Kernel::Plane_3*> > edge_supporting_planes
		= mesh.template property_map<Polygon_mesh::Edge_index, std::set<const Kernel::Plane_3*> >("e:supp_plane").first;//属性映射_每条边所在平面的集合


	Polygon_mesh::Halfedge_index cur = mesh.halfedge(face);//给定面的第一个半边的索引
	Polygon_mesh::Halfedge_index end = cur;//记录 cur ，用于结束循环
	do {
		//获取半边和其所在的平面
		Polygon_mesh::Edge_index ed = mesh.edge(cur);
		const std::set<const Kernel::Plane_3*>& supporting_planes = edge_supporting_planes[ed];
		if (supporting_planes.find(plane_cutting) != supporting_planes.end())
			return;//边的平面集合中包含切割平面，则说明当前边在切割平面上，退出函数

		//获取边的两个顶点及坐标
		Polygon_mesh::Vertex_index s_vd = mesh.source(cur);
		Polygon_mesh::Vertex_index t_vd = mesh.target(cur);
		const Kernel::Point_3& s = mesh.points()[s_vd];
		const Kernel::Point_3& t = mesh.points()[t_vd];

		//判断两个顶点相对于切割平面的位置
		CGAL::Oriented_side s_side = plane_cutting->oriented_side(s);
		CGAL::Oriented_side t_side = plane_cutting->oriented_side(t);
		if (t_side == CGAL::ON_ORIENTED_BOUNDARY)//目标顶点在plane_cutting上
		{
			if (s_side == CGAL::ON_ORIENTED_BOUNDARY)//源顶点在plane_cutting上
				return;//边在plane_cutting上，不需计算交点，退出函数
			else
				existing_vts.push_back(t_vd);//目标点为交点，添加到已有顶点集
		}
		else
			if ((s_side == CGAL::ON_POSITIVE_SIDE && t_side == CGAL::ON_NEGATIVE_SIDE) ||
				(s_side == CGAL::ON_NEGATIVE_SIDE && t_side == CGAL::ON_POSITIVE_SIDE))//源顶点和目标顶点在plane_cutting的两侧
			{
				//当前边与plane_cutting在内部相交
				//计算两顶点到plane_cutting的距离的平方
				FT s_sdist = CGAL::squared_distance(*plane_cutting, s);
				FT t_sdist = CGAL::squared_distance(*plane_cutting, t);

				if (s_sdist <= CGAL::snap_squared_distance_threshold<FT>())//源顶点到plane_cutting的距离小于阈值，说明plane_cutting在源顶点处切割
					existing_vts.push_back(s_vd);//源顶点添加到已有顶点集
				else
					if (t_sdist <= CGAL::snap_squared_distance_threshold<FT>())//目标顶点到plane_cutting的距离小于阈值
						existing_vts.push_back(t_vd);//目标顶点添加到已有顶点集
					else
					{
						const Kernel::Plane_3* plane1 = *(supporting_planes.begin());//边所在平面中的第一个平面
						const Kernel::Plane_3* plane2 = *(supporting_planes.rbegin());//边所在平面中的最后一个平面
						const Kernel::Plane_3* plane3 = const_cast<const Kernel::Plane_3*>(plane_cutting);//获取plane_cutting

						//如果切割平面和当前边所在的两个平面不同，则三个平面可构成三重交点
						if (plane3 != plane1 && plane3 != plane2)
						{
							sort_increasing(plane1, plane2, plane3);//三个平面按指针升序的顺序排序，便于查询
							//const Kernel::Point_3* p = query_intersection(plane1, plane2, plane3);//查询三个平面的交点
							//2023.11.5 修改——————
							//计算三个平面的交集，如果交集的对象可以转换为一个Kernel::Point_3类型的指针，则三个平面相交于一点
							CGAL::Object obj = CGAL::intersection(*plane1, *plane2, *plane3);
							Kernel::Point_3* p;
							if (const Kernel::Point_3* pt = CGAL::object_cast<Kernel::Point_3>(&obj))
							{
								//复制三重交点的坐标，存储点
								p = new Kernel::Point_3(*pt);
								_triplet_intersections[plane1][plane2][plane3] = p;
							}
							else
								p = nullptr;
							//————————————
							if (p)//查询到有交点
							{
								if (CGAL::squared_distance(*p, s) <= CGAL::snap_squared_distance_threshold<FT>())//交点和源顶点的距离小于阈值，交点视为源顶点
									existing_vts.push_back(s_vd);
								else if (CGAL::squared_distance(*p, t) <= CGAL::snap_squared_distance_threshold<FT>())//交点和目标顶点的距离小于阈值，视为目标顶点
									existing_vts.push_back(t_vd);//点存入已有顶点集
								else
									new_vts.push_back(EdgePos(ed, p));//将交点与当前边构成EdgePos对象，加入到新生成的顶点集
							}
							else
								std::cerr << "Fatal error: should have intersection" << std::endl;//错误_不存在交点
						}
						else
							std::cerr << "Fatal error: should not have duplicated planes." << std::endl;//错误_出现重复平面
					}
			}
			else
			{
				// Nothing needs to do here, we will test the next edge
			}


		cur = mesh.next(cur);

	} while (cur != end);

}

bool HypothesisPlaneProcess::halfedge_exists(Polygon_mesh::Vertex_index v1, Polygon_mesh::Vertex_index v2, const Polygon_mesh & mesh)
{
	Polygon_mesh::Halfedge_index h = mesh.halfedge(v1);
	Polygon_mesh::Halfedge_index end = h;
	do {
		Polygon_mesh::Halfedge_index opp = mesh.opposite(h);
		if (mesh.target(opp) == v2)
			return true;
		h = mesh.prev(opp);
	} while (h != end);
	return false;
}

Polygon_mesh::Halfedge_index HypothesisPlaneProcess::split_edge(Polygon_mesh & mesh, const EdgePos & ep, const Kernel::Plane_3 * cutting_plane)
{
	// 这个函数将分割由 'ep' 表示的边
	// - 给新生成的边分配支撑平面
	// - 返回指向新生成顶点的半边索引
	// 函数内部使用了Euler split_edge()操作

	//1、获取mesh的属性映射
	typename Polygon_mesh::template Property_map<Polygon_mesh::Edge_index, std::set<const Kernel::Plane_3*> > edge_supporting_planes
		= mesh.template property_map<Polygon_mesh::Edge_index, std::set<const Kernel::Plane_3*> >("e:supp_plane").first;//属性映射_每个边所在的支撑面

	typename Polygon_mesh::template Property_map<Polygon_mesh::Vertex_index, std::set<const Kernel::Plane_3*> > vertex_supporting_planes
		= mesh.template property_map<Polygon_mesh::Vertex_index, std::set<const Kernel::Plane_3*> >("v:supp_plane").first;//属性映射_每个顶点所在的支撑面

	std::set<const Kernel::Plane_3*> sfs = edge_supporting_planes[ep.edge];//引用_分割边所在的平面//不使用常量引用，因为会在分割后失效
	CGAL_assertion(sfs.size() == 2);//断言_集合仅有两个元素（分割的两个面）


	//2、分割操作
	Polygon_mesh::Halfedge_index h = CGAL::Euler::split_edge(mesh.halfedge(ep.edge), mesh);//在要分割的边上插入一个新顶点，并返回指向新顶点的半边索引
	if (h == Polygon_mesh::null_halfedge())
		return h;//分割失败，返回空的半边索引

	Polygon_mesh::Vertex_index v = mesh.target(h);//新顶点对应的顶点索引
	if (v == Polygon_mesh::null_vertex())
		return Polygon_mesh::null_halfedge();//索引不存在，返回空的半边索引

	typename Polygon_mesh::template Property_map<Polygon_mesh::Vertex_index, Kernel::Point_3>& coords
		= mesh.points();//属性映射_每个顶点坐标
	coords[v] = *ep.pos;//设置新顶点的坐标为EdgePos对象中存储的交点坐标

	//获取新生成的两条边索引和其所在的平面
	Polygon_mesh::Edge_index e1 = mesh.edge(h);
	edge_supporting_planes[e1] = sfs;
	Polygon_mesh::Edge_index e2 = mesh.edge(mesh.next(h));
	edge_supporting_planes[e2] = sfs;

	//分割边所在的平面和切割面加入到新顶点的平面集合中
	vertex_supporting_planes[v] = sfs;
	vertex_supporting_planes[v].insert(cutting_plane);
	CGAL_assertion(vertex_supporting_planes[v].size() == 3);//断言_新顶点的平面集合中有三个元素

	return h;
}

std::vector<Polygon_mesh::Face_index> HypothesisPlaneProcess::split_plane(Polygon_mesh::Face_index face, const Kernel::Plane_3 * cutting_plane, Polygon_mesh & mesh)
{
	std::vector<Polygon_mesh::Face_index> new_faces;//容器_切割后产生的新面的索引

	//1、定义属性映射
	typename Polygon_mesh::template Property_map<Polygon_mesh::Face_index, const Kernel::Plane_3*> face_supporting_planes =
		mesh.template property_map<Polygon_mesh::Face_index, const Kernel::Plane_3*>("f:supp_plane").first;//属性映射_每个面的支撑面

	const Kernel::Plane_3* supporting_plane = face_supporting_planes[face];
	if (supporting_plane == cutting_plane)
		return new_faces;//给定面所在的平面和切割平面相同，则不需要切割，返回空向量


	typename Polygon_mesh::template Property_map<Polygon_mesh::Face_index, CGAL::internal::Planar_segment<Kernel>*> face_supporting_segments =
		mesh.template property_map<Polygon_mesh::Face_index, CGAL::internal::Planar_segment<Kernel>*>("f:supp_segment").first;//属性映射_每个面的切割面

	typename Polygon_mesh::template Property_map<Polygon_mesh::Edge_index, std::set<const Kernel::Plane_3*> > edge_supporting_planes =
		mesh.template property_map<Polygon_mesh::Edge_index, std::set<const Kernel::Plane_3*> >("e:supp_plane").first;//属性映射_每条边的支撑面
	CGAL::internal::Planar_segment<Kernel>* supporting_segment = face_supporting_segments[face];

	//2023.11.23new:添加新的属性映射以表示平面所属的区域
	Polygon_mesh::Property_map<Polygon_mesh::Face_index, std::size_t> face_region_map
		= mesh.property_map<Polygon_mesh::Face_index, std::size_t>("f:region_map").first;//每个面的所属区域
	size_t index_region = face_region_map[face];
	//—————————————————————————


	//2、计算交点，分情况分析切割情况
	std::vector<Polygon_mesh::Vertex_index> existing_vts;//容器_与切割面相交的已有顶点
	std::vector<EdgePos> new_vts;//容器_与切割面相交的新生成顶点
	compute_intersections(mesh, face, cutting_plane, existing_vts, new_vts);//计算给定面和切割平面的交点


	if (existing_vts.size() + new_vts.size() != 2)
		return new_faces;//交点数没有两个，则切割失败，返回空向量
	else if (existing_vts.size() == 2)
	{
		if (existing_vts[0] == existing_vts[1])
			return new_faces;//两交点为同一个顶点，返回空向量

		if (halfedge_exists(existing_vts[0], existing_vts[1], mesh))
			return new_faces;//两交点之间已存在一条边，返回空向量
	}

	//确定两交点的半边
	Polygon_mesh::Halfedge_index h0 = Polygon_mesh::null_halfedge();
	Polygon_mesh::Halfedge_index h1 = Polygon_mesh::null_halfedge();
	if (existing_vts.size() == 2)
	{
		//cutting_plane 在两个已有顶点处切割给定面（而不是一条边）
		h0 = mesh.halfedge(existing_vts[0]);//第一个已有顶点所在半边为h0
		h1 = mesh.halfedge(existing_vts[1]);//第二个已有顶点所在半边为h1
	}
	else if (existing_vts.size() == 1)
	{
		//cutting_plane 在一个已有顶点和一条边的内部切割给定面
		h0 = mesh.halfedge(existing_vts[0]);
		h1 = split_edge(mesh, new_vts[0], cutting_plane);//以新生成顶点的边为h1
	}
	else if (new_vts.size() == 2)
	{
		//cutting_plane 在两条边的内部切割给定面
		h0 = split_edge(mesh, new_vts[0], cutting_plane);
		h1 = split_edge(mesh, new_vts[1], cutting_plane);
	}
	CGAL_assertion(h0 != Polygon_mesh::null_halfedge());//断言_h0、h1不是空的半边索引
	CGAL_assertion(h1 != Polygon_mesh::null_halfedge());

	//为了分割平面，边h0和h1必须附着在同一个面上
	if (mesh.face(h0) != face)
	{
		//h0所在的面不是给定面
		Polygon_mesh::Halfedge_index end = h0;//记录h0的初始值，用于结束循环
		do {
			h0 = mesh.opposite(mesh.next(h0));//将h0沿着给定面的边界顺时针移动一步
			if (mesh.face(h0) == face)
				break;
		} while (h0 != end);
	}
	CGAL_assertion(mesh.face(h0) == face);//断言_h0所在的面是给定面

	if (mesh.face(h1) != face)
	{
		Polygon_mesh::Halfedge_index end = h1;
		do {
			h1 = mesh.opposite(mesh.next(h1));
			if (mesh.face(h1) == face)
				break;
		} while (h1 != end);
	}
	CGAL_assertion(mesh.face(h1) == face);

	//以h0和h1为参数，在给定面上插入一条新边，返回新边对应的一个半边索引
	Polygon_mesh::Halfedge_index h = CGAL::Euler::split_face(h0, h1, mesh);
	if (h == Polygon_mesh::null_halfedge() || mesh.face(h) == Polygon_mesh::null_face())
	{
		std::cerr << "Fatal error. could not split face" << std::endl;//错误_边插入失败/新边没有对应的有效面
		return new_faces;
	}


	//3、平面分割/赋值属性映射
	//在边h的属性映射_支撑面中加入给定面和切割平面
	Polygon_mesh::Edge_index e = mesh.edge(h);
	edge_supporting_planes[e].insert(supporting_plane);
	edge_supporting_planes[e].insert(cutting_plane);
	CGAL_assertion(edge_supporting_planes[e].size() == 2);//断言_边h所在的平面集合中只有两个元素

	// Now the two faces
	Polygon_mesh::Face_index f1 = mesh.face(h);//获取新边所在的一个半边对应的面索引
	face_supporting_segments[f1] = supporting_segment;//将给定面所在的平面段赋值给f1所在的平面段
	face_supporting_planes[f1] = supporting_plane;//将给定面所在的平面赋值给f1所在的平面
	face_region_map[f1] = index_region;//给定面对应的点云区域赋值给f1
	new_faces.push_back(f1);//将f1加入到容器_新面中

	Polygon_mesh::Face_index f2 = mesh.face(mesh.opposite(h));//获取新边所在的另一个半边对应的面索引
	face_supporting_segments[f2] = supporting_segment;
	face_supporting_planes[f2] = supporting_plane;
	new_faces.push_back(f2);
	face_region_map[f2] = index_region;

	return new_faces;
}
