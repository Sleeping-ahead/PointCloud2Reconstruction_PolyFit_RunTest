#pragma once
#include "DataStruct.h"

#define CGAL_USE_SCIP
#include <CGAL/SCIP_mixed_integer_program_traits.h>

using MIP_Solver = CGAL::SCIP_mixed_integer_program_traits<double>;


class ReconstructProcess
{
public:

	//结构体_平面相交半边的点结构，继承自容器_多边形网络半边索引
	struct Intersection : public std::vector<Polygon_mesh::Halfedge_index>
	{
		const Kernel::Point_3* s;//指针_相交半边的起点坐标
		const Kernel::Point_3* t;//指针_相交半边的终点坐标
	};


	//0、参数初始化
	inline void InitParameter(const double& wt_fitting, const double& wt_coverage, const double& wt_complexity);

	//1、从多边形网络中提取邻接关系，获取相邻面之间的交点
	void GetAdjacency(const Polygon_mesh& mesh_candidate, Planes_intersections map_planes_intersection);

	//2、定义约束方程的条件
	void DefineConstraintFactors(Polygon_mesh& mesh_candidate);

	//2.1、定义约束方程的条件_每个能量项做归一化处理
	void DefineConstraintFactors_v2(Polygon_mesh& mesh_candidate);

	//2.2、定义额外条件_面片到顶部距离
	void DefineConstraintFactor_TopDis(const Polygon_mesh& mesh_candidate);
	
	//启用额外的能量偏好项
	inline void UseExtraFactor(const double wt_roof);

	//3、构建约束方程
	void BuildConstraint(const Polygon_mesh& mesh_candidate);

	//4、求解线性约束整数方程
	void SolveConstraint(Polygon_mesh& mesh_candidate, Polygon_mesh & model);

	/*main、使用SCIP求解约束整数程序，重建三维模型表面*/
	bool Reconstruct(Polygon_mesh& mesh_candidate, Polygon_mesh& model,
		const Planes_intersections& map_planes_intersection,
		double wt_fitting = 0.43f, /* = 0.43 ，拟合项*/
		double wt_coverage = 0.27f, /* = 0.27 ，点覆盖项*/
		double wt_complexity = 0.30f); /* = 0.30 ，复杂度项*/

	//5、合并属同一平面的面片，简化显示//TODO:尝试误差为零的二次误差指标的表面简化
	void MergePlane(Polygon_mesh& model);


private:

	std::vector<Intersection> _adjacency;//候选面集的半边集合

	//权重系数
	double _wt_fitting;//拟合项
	double _wt_coverage;//点覆盖项
	double _wt_complexity;//复杂度项

	MIP_Solver solver;//MixedIntegerProgramTraits = MIP_Solver
	std::vector<MIP_Solver::Variable*> variables;//约束条件的变量个数
	
	std::unordered_map<const Intersection*, std::size_t> edge_usage_status;//无序映射_交叉边索引，表示该交叉边被保留或移除
	std::unordered_map<const Intersection*, std::size_t> edge_sharp_status;//无序映射_锐利边索引，表示其是否为锐利边

	//额外的能量偏好项设置
	bool _isUseRoofObj = false;
	double _wt_roof;//顶面偏好项

	//5.1、根据指定的面下标来拟合平面
	Kernel::Plane_3 PlaneFitting(Polygon_mesh mesh, Polygon_mesh::Face_index face);

	//2.2、查看各能量函数值
	void CheckObjectiveValue(const Polygon_mesh& mesh_candidate);


};

inline void ReconstructProcess::InitParameter(const double& wt_fitting, const double& wt_coverage, const double& wt_complexity)
{
	//清空私有变量
	solver.clear();
	edge_usage_status.clear();
	edge_sharp_status.clear();
	variables.clear();

	//权重赋值
	_wt_fitting = wt_fitting;
	_wt_coverage = wt_coverage;
	_wt_complexity = wt_complexity;

}

inline void ReconstructProcess::UseExtraFactor(const double wt_roof)
{
	_isUseRoofObj = true;
	_wt_roof = wt_roof;
}

