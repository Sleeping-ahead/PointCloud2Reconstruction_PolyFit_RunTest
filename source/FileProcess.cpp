#include "FileProcess.h"

std::vector<std::string> FileProcess::GetFileList(const std::string & path)
{
	namespace fs = std::experimental::filesystem;
	std::vector<std::string> files;//文件名列表

	fs::path _path(path);
	if (!fs::exists(_path) || !fs::is_directory(_path))
	{
		std::cout << "错误：此路径不存在 " << path << std::endl;
		return std::vector<std::string>();
	}

	fs::directory_iterator list(_path);	//文件入口容器
	for (const auto& it : list)
	{
		std::string path_file = path + it.path().filename().string();//通过文件入口（it）获取path对象，再得到path对象的文件名，将之输出
		files.push_back(path_file);
	}

	return files;
}

std::string FileProcess::GetFileName(const std::string& path_file)
{
	
	std::size_t lastSlash = path_file.find_last_of("/\\");//字符串中最后一个正or反斜杠的位置索引
														  //适应多系统
	std::size_t lastDot = path_file.find_last_of(".");//最后一个点号的位置索引

	if (lastSlash == std::string::npos) 
		lastSlash = 0;//没有斜杠，从字符串开始处提取
	else
		lastSlash += 1;//从斜杠后提取

	if (lastDot == std::string::npos || lastDot < lastSlash)
		lastDot = path_file.length();//没有点号，或点号在最后一个斜杠之前，设提取结束处在字符串末尾

	return path_file.substr(lastSlash, lastDot - lastSlash);
}

int FileProcess::File2PointCloud(const std::string & path_input, PointCloud & cloud)
{
	std::string suffix = path_input.substr(path_input.find_last_of(".") + 1);//文件名后缀

	if (suffix == "txt" || suffix == "xyz")
	{
		std::ifstream stream(path_input, std::ios_base::binary);//以名文方式打开指定路径的文件
		if (!stream)
		{
			std::cerr << "错误：无法读取到数据所在路径" << path_input << std::endl;
			return 0;
		}
		stream >> cloud;
	}
	else if (suffix == "ply")
	{
		if (!CGAL::IO::read_PLY(path_input, cloud))
		{
			std::cerr << "错误：无法读取到数据所在路径" << std::endl;
			return 0;
		}
	}

	std::cout << "文件已读取，路径：" << path_input << std::endl;
		//<< "输入数据有 " << cloud.size() << " 个点" << std::endl;

	return static_cast<int>(cloud.size());
}

void FileProcess::GiveCloudRandomColor(PointCloud & cloud)
{
	//设置随机颜色
	const unsigned char r = static_cast<unsigned char>(std::rand() % 256);
	const unsigned char g = static_cast<unsigned char>(std::rand() % 256);
	const unsigned char b = static_cast<unsigned char>(std::rand() % 256);

	//对点云中的每个点赋予rgb颜色
	cloud.add_property_map<unsigned char>("red", r).first;
	cloud.add_property_map<unsigned char>("green", g).first;
	cloud.add_property_map<unsigned char>("blue", b).first;

}

void FileProcess::SaveCloud(const std::string & path_output, const std::string & file_name, const std::vector<Kernel::Point_3>& data)
{
	std::string path_entire = path_output + file_name + ".txt";
	std::cout << "保存文件：" << path_entire << std::endl;

	std::ofstream ofs;
	ofs.flags(std::ios::fixed);
	//ofs.precision(10);//保留小数点后3位
	ofs.open(path_entire, std::ios::out);

	//std::srand(std::time(0));//初始化随机数生成器，需设置在函数外
	//生成随机RGB值
	//int r = std::rand() % 256;//随机生成0到255之间的整数
	//int g = std::rand() % 256;
	//int b = std::rand() % 256;

	for (int i = 0; i < data.size(); ++i)
		ofs << data[i].x() << " "
		<< data[i].y() << " "
		<< data[i].z() << " "
		//<< r << " " << g << " " << b << " "
		<< std::endl;

	ofs.close();

}

void FileProcess::SaveCloud(const std::string & path_output, const std::string & file_name, const PointCloud & data)
{
	std::string path_entire = path_output + file_name + ".ply";
	std::cout << "保存文件：" << path_entire << std::endl;

	std::ofstream ofs(path_entire);
	ofs.precision(10);//保留小数点后3位//避免点云复位后丢失精度
	ofs << data;
	ofs.close();
}

bool FileProcess::SaveOBJ(const std::string & path_output, const std::string & file_name, const CGAL::Surface_mesh<Kernel::Point_3>& mesh)
{
	if (mesh.faces().size() == 0)
		return false;

	//获取时间
	//auto now = std::chrono::system_clock::now();//获取当前时间点，使用系统时钟	
	//std::time_t now_time = std::chrono::system_clock::to_time_t(now);//转换为time_t类型，以便使用传统的时间函数
	//std::string timeStr = std::ctime(&now_time);
	
	std::string path_entire = path_output + file_name;
	//输出文件.obj
	std::ofstream out(path_entire + ".obj");//打开一个输出文件流		
	//out << "# file written at " + timeStr + "\n";//文件头
	out << "# " << mesh.number_of_vertices() << " vertices\n"
		<< "# " << mesh.number_of_faces() << " facets\n"
		<< "# " << mesh.number_of_edges() << " halfedges\n\n";//顶点、面、边信息

	out << "# List of geometric vertices\n";//顶点坐标
	for (auto v : vertices(mesh))
		out << "v " << mesh.point(v) << "\n";

	out << "# List of polygonal face elements\n";
	int num_idx = 0;
	for (auto f : faces(mesh))
	{
		out << "f";
		for (auto v : vertices_around_face(mesh.halfedge(f), mesh))
			out << " " << get(CGAL::vertex_index, mesh)[v] + 1; //面的顶点索引。注意.obj文件中的索引从1开始

		out << "\n\n";
	}

	num_idx = 0;
	out.close();// 关闭文件流


	return true;
}

bool FileProcess::SaveColorOBJ(const std::string & path_output, const std::string & file_name, const CGAL::Surface_mesh<Kernel::Point_3>& mesh)
{
	if (mesh.faces().size() == 0)
		return false;

	//获取时间
	//auto now = std::chrono::system_clock::now();//获取当前时间点，使用系统时钟	
	//std::time_t now_time = std::chrono::system_clock::to_time_t(now);//转换为time_t类型，以便使用传统的时间函数
	//std::string timeStr = std::ctime(&now_time);

	std::string path_entire = path_output + file_name;
	std::cout << "保存文件：" << path_entire + ".obj" << std::endl;
	//1、输出文件.obj
	std::ofstream out(path_entire + ".obj");//打开一个输出文件流		
	//out << "# file written at " + timeStr + "\n";//文件头
	out << "# " << mesh.number_of_vertices() << " vertices\n"
		<< "# " << mesh.number_of_faces() << " facets\n"
		<< "# " << mesh.number_of_edges() << " halfedges\n\n";//顶点、面、边信息

	out << "mtllib " + file_name + ".mtl\n\n";//设置引用的.mtl文件

	out << "# List of geometric vertices\n";//顶点坐标
	for (auto v : vertices(mesh))
		out << "v " << mesh.point(v) << "\n";

	out << "# List of polygonal face elements\n";
	int num_idx = 0;
	for (auto f : faces(mesh))
	{
		out << "g group" << num_idx << "\n"
			<< "usemtl color" << num_idx << "\n";//定义分组和材质
		num_idx++;

		out << "f";
		for (auto v : vertices_around_face(mesh.halfedge(f), mesh))
			out << " " << get(CGAL::vertex_index, mesh)[v] + 1; //面的顶点索引。注意.obj文件中的索引从1开始

		out << "\n\n";
	}

	num_idx = 0;
	out.close();// 关闭文件流


	//2、输出文件.mtl
	std::ofstream out_mtl(path_entire + ".mtl");//打开一个输出文件流		
	for (const Polygon_mesh::Face_index& f : faces(mesh))
	{
		out_mtl << "newmtl color" << num_idx << "\n";//定义材质
		num_idx++;

		CGAL::Color c = CGAL::Color(rand() % 256, rand() % 256, rand() % 256); //生成随机颜色
		out_mtl << "Ka " << c.red() / 255.0 << " " << c.green() / 255.0 << " " << c.blue() / 255.0 << "\n";	//环境反射颜色，范围为[0,1]
		out_mtl << "Kd " << c.red() / 255.0 << " " << c.green() / 255.0 << " " << c.blue() / 255.0 << "\n";	//漫反射颜色，范围为[0,1]
		out_mtl << "Ks " << 0.200000 << " " << 0.200000 << " " << 0.200000 << "\n";							//镜面反射颜色，范围为[0,1]
		//out_mtl << "Ks " << c.red() / 255.0 << " " << c.green() / 255.0 << " " << c.blue() / 255.0 << "\n";//耀眼的金属反光

		out_mtl << "\n\n";
	}

	out_mtl.close();// 关闭文件流


	return true;
}

bool FileProcess::SaveColorOFF(const std::string & path_output, const std::string & file_name, const CGAL::Surface_mesh<Kernel::Point_3>& mesh)
{
	if (mesh.faces().size() == 0)
		return false;

	//.off输出文件		
	std::string path_entire = path_output + file_name + ".off";
	std::cout << "保存文件：" << path_entire << std::endl;
	std::ofstream out(path_entire);//打开一个输出文件流		
	out << "OFF\n";//文件头
	out << mesh.number_of_vertices() << " "
		<< mesh.number_of_faces() << " "
		<< mesh.number_of_edges() << "\n\n";//顶点、面、边信息

	for (auto v : vertices(mesh))//顶点坐标
		out << mesh.point(v) << "\n";

	for (auto f : faces(mesh))//面的顶点索引和颜色信息
	{
		out << mesh.degree(f) << " ";//面的顶点数
		for (auto v : vertices_around_face(mesh.halfedge(f), mesh))
			out << " " << get(CGAL::vertex_index, mesh)[v];

		CGAL::Color c = CGAL::Color(rand() % 256, rand() % 256, rand() % 256); //生成随机颜色
		out << " " << c << "\n";//颜色信息
	}

	out.close();// 关闭文件流
	return true;
}


