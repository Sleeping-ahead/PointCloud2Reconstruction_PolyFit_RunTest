#include "ReconstructProcess.h"


void ReconstructProcess::GetAdjacency(const Polygon_mesh & mesh_candidate, Planes_intersections map_planes_intersection)
{
	//1、定义、声明
	typename Polygon_mesh::template Property_map<Polygon_mesh::Vertex_index, std::set<const Kernel::Plane_3*> > vertex_supporting_planes
		= mesh_candidate.template property_map<Polygon_mesh::Vertex_index, std::set<const Kernel::Plane_3*> >("v:supp_plane").first;//属性映射_每个点的支撑面

	// 边由它的两个端点表示
	typedef typename std::unordered_map<const Kernel::Point_3*, std::set<Polygon_mesh::Halfedge_index> >	Edge_map;//无序映射_点指针对应半边索引集合，表示边与它的半边
	typedef typename std::unordered_map<const Kernel::Point_3*, Edge_map >									Face_pool;//无序映射_点指针对应Edge_map，表示一组邻接关系
	Face_pool face_pool;


	//2、向 face_pool 中存入非边界边的所在面及边的两个顶点坐标，每个顶点为三重交点
	for (auto h : mesh_candidate.halfedges())
	{
		Polygon_mesh::Face_index f = mesh_candidate.face(h);//获取每个半边所在的面索引
		if (f == Polygon_mesh::null_face())
			continue;//半边没有对应的面，则为边界半边，跳过

		//获取半边的目标顶点+源顶点、其所在的平面
		Polygon_mesh::Vertex_index sd = mesh_candidate.source(h);
		Polygon_mesh::Vertex_index td = mesh_candidate.target(h);
		const std::set<const Kernel::Plane_3*>& set_s = vertex_supporting_planes[sd];
		const std::set<const Kernel::Plane_3*>& set_t = vertex_supporting_planes[td];
		CGAL_assertion(set_s.size() == 3);//断言_点所在的平面集合有三个元素，即顶点是由三个平面相交而成的
		CGAL_assertion(set_t.size() == 3);

		std::vector<const Kernel::Plane_3*> s_planes(set_s.begin(), set_s.end());//将平面集合转存至容器中，方便调取
		CGAL_assertion(s_planes[0] < s_planes[1]);//断言_平面按序排列
		CGAL_assertion(s_planes[1] < s_planes[2]);
		const Kernel::Point_3* s = map_planes_intersection[s_planes[0]][s_planes[1]][s_planes[2]];//获取交点坐标

		std::vector<const Kernel::Plane_3*> t_planes(set_t.begin(), set_t.end());
		CGAL_assertion(t_planes[0] < t_planes[1]);
		CGAL_assertion(t_planes[1] < t_planes[2]);
		const Kernel::Point_3* t = map_planes_intersection[t_planes[0]][t_planes[1]][t_planes[2]];

		if (s > t)
			std::swap(s, t);//源顶点坐标大于目标顶点，则它们的顺序不一致，交换两者的值
		face_pool[s][t].insert(mesh_candidate.halfedge(f));//源顶点和目标顶点坐标为键，半边所在的面索引为值，插入到 face_pool
	}


	//3、使用 _adjacency 转存 face_pool 中的信息
	for (Face_pool::const_iterator it = face_pool.begin(); it != face_pool.end(); ++it)
	{
		const Kernel::Point_3* s = it->first;//源顶点坐标
		const Edge_map& tmp = it->second;
		typename Edge_map::const_iterator cur = tmp.begin();
		for (; cur != tmp.end(); ++cur)
		{
			const Kernel::Point_3* t = cur->first;//目标顶点坐标
			const std::set<Polygon_mesh::Halfedge_index>& faces = cur->second;//集合_半边索引
			Intersection fan;
			fan.s = s;
			fan.t = t;
			fan.insert(fan.end(), faces.begin(), faces.end());
			_adjacency.push_back(fan);
		}
	}

	///输出测试_adjacency
	//std::cout << "输出测试 adjacency数量：" << _adjacency.size() << std::endl;
	//for (int i = 0; i < _adjacency.size(); ++i)
	//{
	//	std::cout << _adjacency[i].s->x() << " " << _adjacency[i].s->y() << " " << _adjacency[i].s->z() << std::endl;
	//	std::cout << _adjacency[i].t->x() << " " << _adjacency[i].t->y() << " " << _adjacency[i].t->z() << std::endl;
	//	std::cout << std::endl;
	//}

}

void ReconstructProcess::DefineConstraintFactors(Polygon_mesh& mesh_candidate)
{
	//1、参数定义/起别名
	//mesh_candidate属性映射
	typename Polygon_mesh::template Property_map<Polygon_mesh::Face_index, std::size_t> face_num_supporting_points =
		mesh_candidate.template property_map<Polygon_mesh::Face_index, std::size_t>("f:num_supporting_points").first;//属性映射_面的支撑点数量

	typename Polygon_mesh::template Property_map<Polygon_mesh::Face_index, FT> face_areas =
		mesh_candidate.template property_map<Polygon_mesh::Face_index, FT>("f:face_area").first;//属性映射_面的面积

	typename Polygon_mesh::template Property_map<Polygon_mesh::Face_index, FT> face_covered_areas =
		mesh_candidate.template property_map<Polygon_mesh::Face_index, FT>("f:covered_area").first;//属性映射_面的点覆盖面积

	typename Polygon_mesh::template Property_map<Polygon_mesh::Face_index, std::size_t> face_indices =
		mesh_candidate.template add_property_map<Polygon_mesh::Face_index, std::size_t>("f:index").first;//属性映射_面索引

	//支撑点总数
	double total_points = 0.0;
	std::size_t idx = 0;
	for (auto f : mesh_candidate.faces())
	{
		total_points += face_num_supporting_points[f];
		face_indices[f] = idx;//给每个面附上新的索引值
		++idx;
	}

	std::size_t num_faces = mesh_candidate.number_of_faces();//候选表面的面数量
	std::size_t num_edges(0);//交叉边数量，初始值为0

	//二进制变量定义：
	//x[0] ... x[num_faces - 1]:输入面的二进制标签
	//x[num_faces] ... x[num_faces + num_edges - 1]:相交边的二进制标签（保留或不保留）
	//x[num_faces + num_edges] ...  x[num_faces + num_edges + num_edges]:角边的二进制标签（非锐边）


	//2、计算约束条件的值
	// edge_usage_status 赋值
	for (std::size_t i = 0; i < _adjacency.size(); ++i)
	{
		//遍历所有含邻接关系的边，并将包含4个半边的交叉边索引添加到 edge_usage_status
		//该边由4个候选面相交而成
		const Intersection& fan = _adjacency[i];
		if (fan.size() == 4)
		{
			//交叉边索引 为键，候选面数量+交叉边数量 为值
			std::size_t var_idx = num_faces + num_edges;
			edge_usage_status[&fan] = var_idx;
			++num_edges;
		}
	}

	// variables 定义、赋值
	std::size_t total_variables = num_faces + num_edges + num_edges;//总变量数量 = 候选面数量 + 2 * 交叉边数量
	variables = solver.create_variables(total_variables);//容器_变量个数为 total_variables
	for (std::size_t i = 0; i < total_variables; ++i)
	{
		//为每个变量设置类型为 BINARY二进制，只能取0或1
		MIP_Solver::Variable* v = variables[i];
		v->set_variable_type(MIP_Solver::Variable::BINARY);
	}


	//3、添加MIP_Solver的 objective（各个能量项？）
	//候选面集的顶点坐标及其包围盒参数
	std::vector<Kernel::Point_3> vertices(mesh_candidate.number_of_vertices());
	idx = 0;
	for (auto v : mesh_candidate.vertices())
	{
		vertices[idx] = mesh_candidate.points()[v];
		++idx;
	}
	const Kernel::Iso_cuboid_3& box = CGAL::bounding_box(vertices.begin(), vertices.end());//计算 vertices 点集的最小包围盒
	FT dx = box.xmax() - box.xmin();//包围盒xyz方向上的长度
	FT dy = box.ymax() - box.ymin();
	FT dz = box.zmax() - box.zmin();
	FT box_area = FT(2.0) * (dx * dy + dy * dz + dz * dx);//包围盒表面积

	//选择更好的尺度
	double coeff_data_fitting = _wt_fitting;//拟合项的权重系数
	double coeff_coverage = total_points * _wt_coverage / box_area;//点覆盖项的权重系数
	double coeff_complexity = total_points * _wt_complexity / double(_adjacency.size());//复杂度项的权重系数

	///输出测试
	//std::cout << coeff_data_fitting << "  " << coeff_coverage << "  " << coeff_complexity << std::endl;

	MIP_Solver::Linear_objective * objective = solver.create_objective(MIP_Solver::Linear_objective::MINIMIZE);//指针_线性目标函数对象，并设置类型为最小化

	// edge_sharp_status 赋值
	std::size_t num_sharp_edges = 0;
	for (std::size_t i = 0; i < _adjacency.size(); ++i)
	{
		const Intersection& fan = _adjacency[i];
		if (fan.size() == 4)
		{
			std::size_t var_idx = num_faces + num_edges + num_sharp_edges;//尖锐边数量可能小于交叉边，部分variables全程0值？
			edge_sharp_status[&fan] = var_idx;

			//累计复杂度项
			objective->add_coefficient(variables[var_idx], coeff_complexity);//为每个MIP变量添加系数
																			 //系数 = 复杂度项的权重，表示复杂度项是最小化目标（联系157行来理解）
			++num_sharp_edges;
		}
	}
	CGAL_assertion(num_edges == num_sharp_edges);

	for (auto f : mesh_candidate.faces())
	{
		std::size_t var_idx = face_indices[f];

		//累计数据拟合项
		std::size_t num = face_num_supporting_points[f];//面的支撑点数量
		objective->add_coefficient(variables[var_idx], -coeff_data_fitting * num);//为每个MIP变量添加系数
																				  //系数=拟合能量项的权重*（-支撑点数），表示拟合能量项是最小化目标
		//累计点覆盖项
		double uncovered_area = (face_areas[f] - face_covered_areas[f]);//面的未覆盖面积
		objective->add_coefficient(variables[var_idx], coeff_coverage * uncovered_area);
	}

	//添加约束：与边关联的面数必须为2或0
	std::size_t var_edge_used_idx = 0;//记录当前交叉边对应的变量索引
	for (std::size_t i = 0; i < _adjacency.size(); ++i)
	{
		MIP_Solver::Linear_constraint* c = solver.create_constraint(0.0, 0.0);//线性约束_约束上界和下界 = 0，表示无约束条件约束(无取值范围限制
		const Intersection& fan = _adjacency[i];
		for (std::size_t j = 0; j < fan.size(); ++j)
		{
			//对面添加约束
			Polygon_mesh::Face_index f = mesh_candidate.face(fan[j]);
			std::size_t var_idx = face_indices[f];
			c->add_coefficient(variables[var_idx], 1.0);//表示约束中包含变量 variables[var_idx]
		}

		//交叉边由四个候选面相交而成
		if (fan.size() == 4)
		{
			//对交叉边添加约束
			std::size_t var_idx = num_faces + var_edge_used_idx;
			c->add_coefficient(variables[var_idx], -2.0);//表示约束中包含这个变量，并且它与其他变量之和必须等于零？
			++var_edge_used_idx;
		}
		else
		{
			//交叉边为边界边，不添加任何系数，因为生成的模型不为开放的表面模型
		}
	}

}

void ReconstructProcess::DefineConstraintFactors_v2(Polygon_mesh & mesh_candidate)
{
	//1、参数定义/起别名
	//mesh_candidate属性映射
	typename Polygon_mesh::template Property_map<Polygon_mesh::Face_index, std::size_t> face_num_supporting_points =
		mesh_candidate.template property_map<Polygon_mesh::Face_index, std::size_t>("f:num_supporting_points").first;//属性映射_面的支撑点数量

	typename Polygon_mesh::template Property_map<Polygon_mesh::Face_index, FT> face_areas =
		mesh_candidate.template property_map<Polygon_mesh::Face_index, FT>("f:face_area").first;//属性映射_面的面积

	typename Polygon_mesh::template Property_map<Polygon_mesh::Face_index, FT> face_covered_areas =
		mesh_candidate.template property_map<Polygon_mesh::Face_index, FT>("f:covered_area").first;//属性映射_面的点覆盖面积

	typename Polygon_mesh::template Property_map<Polygon_mesh::Face_index, std::size_t> face_indices =
		mesh_candidate.template add_property_map<Polygon_mesh::Face_index, std::size_t>("f:index").first;//属性映射_面索引

	//支撑点总数
	double total_points = 0.0;
	std::size_t idx = 0;
	for (auto f : mesh_candidate.faces())
	{
		total_points += face_num_supporting_points[f];
		face_indices[f] = idx;//给每个面附上新的索引值
		++idx;
	}

	std::size_t num_faces = mesh_candidate.number_of_faces();//候选表面的面数量
	std::size_t num_edges(0);//交叉边数量，初始值为0

	//二进制变量定义：
	//x[0] ... x[num_faces - 1]:输入面的二进制标签
	//x[num_faces]…x[num_faces] ... x[num_faces + num_edges - 1]:相交边的二进制标签（保留或不保留）
	//x[num_faces + num_edges] ...  x[num_faces + num_edges + num_edges]:角边的二进制标签（非锐边）


	//2、计算约束条件的值
	// edge_usage_status 赋值
	for (std::size_t i = 0; i < _adjacency.size(); ++i)
	{
		//遍历所有含邻接关系的边，并将包含4个半边的交叉边索引添加到 edge_usage_status
		//该边由4个候选面相交而成
		const Intersection& fan = _adjacency[i];
		if (fan.size() == 4)
		{
			//交叉边索引 为键，候选面数量+交叉边数量 为值
			std::size_t var_idx = num_faces + num_edges;
			edge_usage_status[&fan] = var_idx;
			++num_edges;
		}
	}

	// variables 定义、赋值
	std::size_t total_variables = num_faces + num_edges + num_edges;//总变量数量 = 候选面数量 + 2 * 交叉边数量
	variables = solver.create_variables(total_variables);//容器_变量个数为 total_variables
	for (std::size_t i = 0; i < total_variables; ++i)
	{
		//为每个变量设置类型为 BINARY二进制，只能取0或1
		MIP_Solver::Variable* v = variables[i];
		v->set_variable_type(MIP_Solver::Variable::BINARY);
	}


	//3、添加MIP_Solver的 objective（各个能量项？）
	//候选面集的顶点坐标及其包围盒参数
	std::vector<Kernel::Point_3> vertices(mesh_candidate.number_of_vertices());
	idx = 0;
	for (auto v : mesh_candidate.vertices())
	{
		vertices[idx] = mesh_candidate.points()[v];
		++idx;
	}
	const Kernel::Iso_cuboid_3& box = CGAL::bounding_box(vertices.begin(), vertices.end());//计算 vertices 点集的最小包围盒
	FT dx = box.xmax() - box.xmin();//包围盒xyz方向上的长度
	FT dy = box.ymax() - box.ymin();
	FT dz = box.zmax() - box.zmin();
	FT box_area = FT(2.0) * (dx * dy + dy * dz + dz * dx);//包围盒表面积


	MIP_Solver::Linear_objective * objective = solver.create_objective(MIP_Solver::Linear_objective::MINIMIZE);//指针_线性目标函数对象，并设置类型为最小化

	// edge_sharp_status 赋值
	std::size_t num_sharp_edges = 0;
	for (std::size_t i = 0; i < _adjacency.size(); ++i)
	{
		const Intersection& fan = _adjacency[i];
		if (fan.size() == 4)//该边为锐边
		{
			std::size_t var_idx = num_faces + num_edges + num_sharp_edges;//该var_idx的计数 ?= 以face_indices[f]获取的计数
			edge_sharp_status[&fan] = var_idx;

			//累计复杂度项
			objective->add_coefficient(variables[var_idx], -_wt_complexity);//为每个MIP变量添加系数
																			//系数 = -复杂度项的权重，权重越高越需要该锐边
			++num_sharp_edges;
		}
	}
	CGAL_assertion(num_edges == num_sharp_edges);


	int num_max = -1;//面的支撑点数量最大值
	for (const Polygon_mesh::Face_index& f : mesh_candidate.faces())
	{
		int num = static_cast<int>(face_num_supporting_points[f]);//面的支撑点数量
		if (num > num_max)
			num_max = num;
	}
	///输出测试_最大的支撑点数
	//std::cout << num_max << std::endl;

	for (auto f : mesh_candidate.faces())
	{
		std::size_t var_idx = face_indices[f];

		//累计数据拟合项
		//std::size_t num = face_num_supporting_points[f];//面的支撑点数量
		double num = static_cast<double>(face_num_supporting_points[f]);//面的支撑点数量
		objective->add_coefficient(variables[var_idx], -_wt_fitting * num / num_max);//为每个MIP变量添加系数
																					 //系数=拟合能量项的权重*（-支撑点数），点数越多越靠近最小化
		///输出测试
		//std::cout << -_wt_fitting * num / num_max << "   "
		//	<< objective->get_coefficient(variables[var_idx]) << std::endl << std::endl;

		//累计点覆盖项
		//double uncovered_area = (face_areas[f] - face_covered_areas[f]);//面的未覆盖面积
		//objective->add_coefficient(variables[var_idx], -_wt_coverage * uncovered_area / face_areas[f]);
		objective->add_coefficient(variables[var_idx], -_wt_coverage * face_covered_areas[f] / face_areas[f]);
	}

	//添加约束：与边关联的面数必须为2或0
	std::size_t var_edge_used_idx = 0;//记录当前交叉边对应的变量索引
	for (std::size_t i = 0; i < _adjacency.size(); ++i)
	{
		MIP_Solver::Linear_constraint* c = solver.create_constraint(0.0, 0.0);//线性约束_约束上界和下界 = 0，表示等式约束
		const Intersection& fan = _adjacency[i];
		for (std::size_t j = 0; j < fan.size(); ++j)
		{
			//对面添加约束
			Polygon_mesh::Face_index f = mesh_candidate.face(fan[j]);
			std::size_t var_idx = face_indices[f];
			c->add_coefficient(variables[var_idx], 1.0);//表示约束中包含变量 variables[var_idx]
		}

		//交叉边由四个候选面相交而成（下标是否存在合理性？）
		if (fan.size() == 4)
		{
			//对交叉边添加约束
			std::size_t var_idx = num_faces + var_edge_used_idx;
			c->add_coefficient(variables[var_idx], -2.0);//表示约束中包含这个变量，并且它与其他变量之和必须等于零
			++var_edge_used_idx;
		}
		else
		{
			//交叉边为边界边，不添加任何系数，因为生成的模型不为开放的表面模型
		}
	}

}

void ReconstructProcess::DefineConstraintFactor_TopDis(const Polygon_mesh & mesh_candidate)
{
	typename Polygon_mesh::template Property_map<Polygon_mesh::Face_index, double> face_roof =
		mesh_candidate.template property_map<Polygon_mesh::Face_index, double>("f:dist_boxTop").first;//属性映射_面片属顶面

	typename Polygon_mesh::template Property_map<Polygon_mesh::Face_index, std::size_t> face_indices =
		mesh_candidate.template property_map<Polygon_mesh::Face_index, std::size_t>("f:index").first;//属性映射_面索引


	//double total_points = 0.0;//支撑点总数
	//std::size_t num_faces = mesh_candidate.number_of_faces();//候选表面的面数量

	MIP_Solver::Linear_objective * objective = solver.objective();

	//对每个能量函数，累计顶面偏好值
	for (auto f : mesh_candidate.faces())
	{
		std::size_t var_idx = face_indices[f];

		double dist_p2boxTop = face_roof[f];
		objective->add_coefficient(variables[var_idx], _wt_roof * dist_p2boxTop);//离包围盒顶面越远，值越大
	}
}

void ReconstructProcess::BuildConstraint(const Polygon_mesh & mesh_candidate)
{
	//为尖锐的边缘添加约束。提出这个限制的解释可以在这里找到:
	// https://user-images.githubusercontent.com/15526536/30185644-12085a9c-942b-11e7-831d-290dd2a4d50c.png


	// mesh_candidate 属性映射
	typename Polygon_mesh::template Property_map<Polygon_mesh::Face_index, const Kernel::Plane_3*> face_supporting_planes =
		mesh_candidate.template property_map<Polygon_mesh::Face_index, const Kernel::Plane_3*>("f:supp_plane").first;//属性映射_面的支撑面

	typename Polygon_mesh::template Property_map<Polygon_mesh::Face_index, std::size_t> face_indices =
		mesh_candidate.template property_map<Polygon_mesh::Face_index, std::size_t>("f:index").first;//属性映射_面索引


	double M = 1.0;//线性规划中的松弛变量
	for (std::size_t i = 0; i < _adjacency.size(); ++i)
	{
		const Intersection& fan = _adjacency[i];//获取容器_邻接关系中的交叉边
		if (fan.size() != 4)
			continue;

		//构建约束：X[var_edge_usage_idx] >= X[var_edge_sharp_idx]。表示如果一个交叉边是锐利的，那么它必须被选择
		MIP_Solver::Linear_constraint* c = solver.create_constraint(0.0);//线性约束对象_下界，表示约束是大于等于0的不等式约束
		std::size_t var_edge_usage_idx = edge_usage_status[&fan];
		c->add_coefficient(variables[var_edge_usage_idx], 1.0);//对使用状态变量添加系数，系数值为1.0
		std::size_t var_edge_sharp_idx = edge_sharp_status[&fan];
		c->add_coefficient(variables[var_edge_sharp_idx], -1.0);//对锐利状态变量添加系数，系数值为-1.0

		for (std::size_t j = 0; j < fan.size(); ++j)
		{
			//获取相邻次序的两个交叉边及其支撑面
			Polygon_mesh::Face_index f1 = mesh_candidate.face(fan[j]);
			const Kernel::Plane_3* plane1 = face_supporting_planes[f1];
			std::size_t fid1 = face_indices[f1];
			for (std::size_t k = j + 1; k < fan.size(); ++k)
			{
				Polygon_mesh::Face_index f2 = mesh_candidate.face(fan[k]);
				const Kernel::Plane_3* plane2 = face_supporting_planes[f2];
				std::size_t fid2 = face_indices[f2];

				if (plane1 != plane2)
				{
					//两个支撑面不共面的情况下，构成约束：X[var_edge_sharp_idx] + M * (3 - (X[fid1] + X[fid2] + X[var_edge_usage_idx])) >= 1，
					//等价于：X[var_edge_sharp_idx] - M * X[fid1] - M * X[fid2] - M * X[var_edge_usage_idx] >= 1 - 3M
					//表示如果一个交叉边是锐利的，那么它所在的四个面中至少有两个被选择，否则就会导致表面不连续
					c = solver.create_constraint(1.0 - 3.0 * M);//线性约束对象_下界，表示约束是大于等于这个值的不等式约束
					c->add_coefficient(variables[var_edge_sharp_idx], 1.0);//线性约束对象_对锐利状态变量添加一个系数，设为1.0
					c->add_coefficient(variables[fid1], -M);//线性约束对象_设置第一面的变量系数为-M
					c->add_coefficient(variables[fid2], -M);//线性约束对象_设置第二面的变量系数为-M
					c->add_coefficient(variables[var_edge_usage_idx], -M);
				}
			}
		}
	}

	///输出测试
	//std::cout << "最后的各个变量项" << std::endl;
	//CheckObjectiveValue(mesh_candidate);

}

void ReconstructProcess::SolveConstraint(Polygon_mesh & mesh_candidate, Polygon_mesh & model)
{
	if (solver.solve())//求解线性规划问题成功
	{
		//1、根据solve中的解，删除 mesh_candidate 中相应的边
		typename Polygon_mesh::template Property_map<Polygon_mesh::Face_index, std::size_t> face_indices =
			mesh_candidate.template property_map<Polygon_mesh::Face_index, std::size_t>("f:index").first;//属性映射_面索引

		const std::vector<double>& X = solver.solution();//线性规划问题的解

		std::vector<Polygon_mesh::Face_index> to_delete;//容器_面索引_记录需要删除的面
		std::size_t f_idx(0);
		for (auto f : mesh_candidate.faces())
		{
			///输出测试2.9
			//std::cout << X[f_idx] << "  ";

			//解的值四舍五入后等于0，说明这个面没有被选择，将该面的索引存入 to_delete
			if (static_cast<int>(std::round(X[f_idx])) == 0)
				to_delete.push_back(f);
			++f_idx;
		}

		//遍历需要删除的面，使用 CGAL::Euler::remove_face() ，根据半边索引删除相应面
		for (std::size_t i = 0; i < to_delete.size(); ++i)
		{
			Polygon_mesh::Face_index f = to_delete[i];
			Polygon_mesh::Halfedge_index h = mesh_candidate.halfedge(f);
			CGAL::Euler::remove_face(h, mesh_candidate);
		}


		//2、设置 mesh_candidate 中的锐利边属性
		typename Polygon_mesh::template Property_map<Polygon_mesh::Edge_index, bool> edge_is_sharp =
			mesh_candidate.template add_property_map<Polygon_mesh::Edge_index, bool>("e:sharp_edges").first;//属性映射_是否锐利边
		for (auto e : mesh_candidate.edges())
			edge_is_sharp[e] = false;//该属性映射默认没有锐利边

		for (std::size_t i = 0; i < _adjacency.size(); ++i)
		{
			const Intersection& fan = _adjacency[i];
			if (fan.size() != 4)
				continue;//跳过边界边

			//获取交叉边的锐利状态变量索引，若该值为1（该边是锐利的）
			std::size_t idx_sharp_var = edge_sharp_status[&fan];
			if (static_cast<int>(X[idx_sharp_var]) == 1)
			{
				for (std::size_t j = 0; j < fan.size(); ++j)
				{
					//获取每个半边的索引和相应的面索引，若面索引非空（该面未被删除）
					Polygon_mesh::Halfedge_index h = fan[j];
					Polygon_mesh::Face_index f = mesh_candidate.face(h);
					if (f != Polygon_mesh::null_face())
					{
						//获取该面，若果解的fid元素四舍五入后值为1
						std::size_t fid = face_indices[f];
						if (static_cast<int>(std::round(X[fid])) == 1)
						{
							Polygon_mesh::Edge_index e = mesh_candidate.edge(h);
							edge_is_sharp[e] = true;//该面被选择，设置该边的 属性映射_锐利 值为真
							break;
						}
					}
				}
			}
		}


		//3、将mesh_candidate数据存入model
		model.clear();
		CGAL::copy_face_graph(mesh_candidate, model);
	}
	else
		std::cerr << "solving the binary program failed" << std::endl;

}

bool ReconstructProcess::Reconstruct(Polygon_mesh & mesh_candidate, Polygon_mesh & model,
	const Planes_intersections & map_planes_intersection,
	double wt_fitting, double wt_coverage, double wt_complexity)
{
	CGAL::Timer timer;//时间计数的实例
	timer.start();

	//0、检验数据有效性
	if (mesh_candidate.num_faces() < 4)
	{
		std::cout << "错误：至少需要4个候选表面才能重建出紧密模型，当前的表面集仅有 " + std::to_string(mesh_candidate.num_faces()) + " 个面";
		return false;
	}
	//参数初始化
	InitParameter(wt_fitting, wt_coverage, wt_complexity);


	//1、计算表面的邻接信息：mesh_candidate的三重交点及其边的所在面
	GetAdjacency(mesh_candidate, map_planes_intersection);


	//2、设置约束条件
	if (_isUseRoofObj)
	{
		DefineConstraintFactors_v2(mesh_candidate);//各能量项归一化处理
		DefineConstraintFactor_TopDis(mesh_candidate);//添加顶面偏好项
	}
	else
		DefineConstraintFactors(mesh_candidate);//原PolyFit方法
	///输出测试_能量项累计值
	//CheckObjectiveValue(mesh_candidate);


	//3、构建约束方程
	BuildConstraint(mesh_candidate);


	//4、求解线性约束整数方程
	SolveConstraint(mesh_candidate, model);

	timer.stop();
	std::cout << "重建模型的面数为 " << model.faces().size() << " 个，运行耗时 " << timer.time() << " 秒" << std::endl;

	return true;
}

void ReconstructProcess::MergePlane(Polygon_mesh & model)
{
	CGAL::Timer timer;//时间计数的实例
	timer.start();

	FT theta = cos(static_cast<FT>(CGAL_PI * 15.0 / FT(180.0)));//阈值_向量的内积（余弦值），大于该值则两平面近似平行

	bool merged = false;//是否有平面被合并
	do
	{
		merged = false;

		//1、遍历model中的边，选取两个邻近的面进行判断
		std::vector<Polygon_mesh::Halfedge_index> index_halfedges;//容器_model的半边集下标
		for (Polygon_mesh::Halfedge_index e : model.halfedges())
			index_halfedges.push_back(e);

		for (int i = 0; i < index_halfedges.size(); ++i)
		{
			if (CGAL::is_border_edge(index_halfedges[i], model))
				continue;//边为多边形上的边界边，跳过

			//获取model的一个半边及其相应的面
			Polygon_mesh::Halfedge_index h1 = index_halfedges[i];
			Polygon_mesh::Face_index face1 = model.face(h1);
			//获取对偶半边和相应的面
			Polygon_mesh::Halfedge_index h2 = model.opposite(h1);
			Polygon_mesh::Face_index face2 = model.face(h2);

			if (face1 == face2)
				continue;//相邻面为同一个面，跳过

			if ((!model.is_valid(face1)) || (!model.is_valid(face2)))
				continue;//相邻的任一个面不存在，跳过


			//2、计算两个平面的正交向量并归一化
			//拟合平面
			Kernel::Plane_3 plane1 = PlaneFitting(model, face1);
			Kernel::Plane_3 plane2 = PlaneFitting(model, face2);
			//计算正交向量，归一化
			Kernel::Vector_3 n1 = plane1.orthogonal_vector();
			Kernel::Vector_3 n2 = plane2.orthogonal_vector();
			CGAL::internal::normalize<FT, Kernel::Vector_3>(n1);
			CGAL::internal::normalize<FT, Kernel::Vector_3>(n2);

			///输出测试_正交向量值
			//std::cout << "这两个平面  " << face1 << "  " << face2 << std::endl;
			//std::cout << plane1 << "  n1 " << n1.x() << " " << n1.y() << " " << n1.z() << " " << std::endl;
			//std::cout << plane2 << "  n2 " << n2.x() << " " << n2.y() << " " << n2.z() << " " << std::endl;
			//std::cout << CGAL::angle(n1, n2) << "  " << n1 * n2 << std::endl << std::endl;

			//如果两个平面正交向量之间的内积大于阈值，合并平面
			if (std::abs(n1 * n2) > std::cos(theta))
			{
				//std::cout << "这两个平面需要合并  " << face1 << "  " << face2 << std::endl;
				CGAL::Euler::join_face(h1, model);
				merged = true;
				break;
			}
		}

	} while (merged);//merged 若为false，则此时所有平面都不适合被合并，跳出循环

	timer.stop();
	std::cout << "合并简化模型的面数为 " << model.faces().size() << " 个，运行耗时 " << timer.time() << " 秒" << std::endl;

}

Kernel::Plane_3 ReconstructProcess::PlaneFitting(Polygon_mesh mesh, Polygon_mesh::Face_index face)
{
	//1、通过遍历 face 的半边，获取 face 的所有顶点
	Polygon_mesh::Halfedge_index edge_cur = mesh.halfedge(face);
	std::vector<Kernel::Point_3> points;
	do
	{
		//获取半边源点的坐标
		Polygon_mesh::Vertex_index v = mesh.source(edge_cur);
		Kernel::Point_3 p = mesh.point(v);
		points.push_back(p);

		edge_cur = mesh.next(edge_cur);
	} while (edge_cur != mesh.halfedge(face));


	//2、创建平面对象，根据点集生成平面
	Kernel::Plane_3 plane;
	CGAL::linear_least_squares_fitting_3(points.begin(), points.end(), plane, CGAL::Dimension_tag<0>());

	///输出测试_点集
	//std::cout << face << std::endl;
	//for (int i = 0; i < points.size(); ++i)
	//	std::cout << points[i] << std::endl;
	//std::cout << std::endl;


	return plane;
}

void ReconstructProcess::CheckObjectiveValue(const Polygon_mesh & mesh_candidate)
{	
	typename Polygon_mesh::template Property_map<Polygon_mesh::Face_index, std::size_t> face_indices =
		mesh_candidate.template property_map<Polygon_mesh::Face_index, std::size_t>("f:index").first;//属性映射_面索引

	MIP_Solver::Linear_objective* objective = solver.objective();

	//获取变量个数
	std::cout << "变量数： " << solver.variables().size() << std::endl;

	//仅获取面片的变量累计值
	for (const Polygon_mesh::Face_index& f : mesh_candidate.faces())
	{
		std::size_t var_idx = face_indices[f];
		std::cout << objective->get_coefficient(variables[var_idx]) << "  ";
	}
	///输出测试_获取各变量累计值
	//for (MIP_Solver::Variable* var : variables)
	//	std::cout << objective->get_coefficient(var) << "  ";
}


//for (size_t var_idx = 0; var_idx < solver.variables().size(); ++var_idx)
//	std::cout << objective->get_coefficient(variables[var_idx]) << "  ";

///输出测试
//std::cout << "三个能量项" << std::endl;
//CheckObjectiveValue(mesh_candidate);

//variables.size() = solver.variables().size()
//&variables[var_idx] = objective->get_coefficient(variables[var_idx])