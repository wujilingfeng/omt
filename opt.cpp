#include "opt.h"
//要保证原测度之和等于现有测度，误差要非常小（关键）
igl::viewer::Viewer viewer;
TMyMesh mesh, mesh1;
Omt myomt(mesh, viewer);

int main()
{
	igl::readOFF("e:/off/test.off", myomt.V, myomt.F);
	Eigen::MatrixXi F,F1;
	Eigen::MatrixXd V,V1;
	Eigen::MatrixXd N_vertices;
	igl::readOFF("e:/off/face.off", V, F);
	igl::per_vertex_normals(V, F, N_vertices);
	OpenMesh::IO::Options opt;
	if (!OpenMesh::IO::read_mesh(mesh, "e:/off/test.off", opt))
	{
		printf_s("读取失败\r\n");
		return 0;
	}
	
	if (!OpenMesh::IO::read_mesh(mesh1, "e:/off/face.off", opt))
	{
		printf_s("读取失败\r\n");
		return 0;
	}
	printf_s("mesh.n-faces %d\r\n", mesh.n_vertices());
	myomt.init();
	printf_s("mesh.n-faces %d\r\n", mesh1.n_vertices());
	HalfedgeHandle he; int i = 0;
	for (PMyMesh::VertexIter vter = mesh1.vertices_begin(); vter != mesh1.vertices_end(); vter++)
	{
		myomt.W.coeffRef(i) = acos(-1) / double(mesh1.n_vertices());
		/*double  sum = 0;
		for (PMyMesh::VertexFaceIter vfter = mesh1.vf_iter(*vter); vfter.is_valid(); vfter++)
		{
			he = mesh.halfedge_handle(*vfter);
			sum += mesh.calc_sector_area(he);
			
		}
		myomt.W.coeffRef(i) = (sum / 3.0);*/
		//mesh.data(*vter).mesure = (sum / 3.0);
		i++;
	}
	
	myomt.W *= (acos(-1) / myomt.W.sum());
	printf_s("W sum %.12lf\r\n",myomt.W.sum());
	
	
	
	myomt.iteration2();
	
	
	myomt.drawedge();

	//myomt.drawaxis();
	
	//viewer.data.set_mesh(myomt.V,myomt.F);
	//viewer.data.set_normals(N_vertices);
	viewer.launch();
	
	printf_s("libo\r\n");
	int n;
	scanf_s("%d", &n);
	return 0;


}
