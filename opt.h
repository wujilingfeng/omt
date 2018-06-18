#pragma once
#include<iostream>
#include <list>
#include<nanogui\screen.h>
#include<igl\readOFF.h>

#include<igl\viewer\Viewer.h>
#include<igl\jet.h>
#include<igl/per_vertex_normals.h>
//#include<nanogui\formhelper.h>
//#include<Eigen\Core>
//#include <Eigen\Eigenvalues>
//#include<OpenMesh\Core\IO\IOManager.hh>
#include<OpenMesh\Core\IO\MeshIO.hh>
#include<OpenMesh\Core\Mesh\PolyMesh_ArrayKernelT.hh>
#include<OpenMesh\Core\Mesh\TriMesh_ArrayKernelT.hh>
#include<OpenMesh\Tools\Dualizer\meshDualT.hh>
#include<Eigen\SparseQR>
//opemesh的遍历默认都是顺时针
using namespace std;
using namespace OpenMesh;
struct MyTraits :OpenMesh::DefaultTraits
{
	VertexTraits{
	
	};
	EdgeTraits{
		

	};
	HalfedgeTraits{

	};
	FaceTraits
	{
	
	};

};
typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> PMyMesh;
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> TMyMesh;
OpenMesh::HalfedgeHandle halfedgeh(VertexHandle v1, VertexHandle v2, PMyMesh &mesh)
{//返回v1起点v2终点的半边
	HalfedgeHandle he = mesh.halfedge_handle(v1), he1 = he;
	//mesh.edge
	if (mesh.from_vertex_handle(he1) != v1)
	{
		cout << "错误：：" << endl;
	}
	do
	{
		he1 = mesh.ccw_rotated_halfedge_handle(he1);
		//he1 = (he1);
	} while (he1 != he && v2 != mesh.to_vertex_handle(he1));
	return he1;
}
PMyMesh *set_meshdual(TMyMesh &primal)
{
	PMyMesh *dual = new PMyMesh();
	FPropHandleT<PMyMesh::VertexHandle> primalToDual;
	EPropHandleT<PMyMesh::VertexHandle> edgev;
	primal.add_property(primalToDual);
	primal.add_property(edgev);
	for (TMyMesh::FaceIter fter = primal.faces_begin(); fter != primal.faces_end(); fter++)
	{
		PMyMesh::Point p(0, 0, 0);
		primal.property(primalToDual, *fter) = dual->add_vertex(p);
	}
	for (TMyMesh::EdgeIter eter = primal.edges_begin(); eter != primal.edges_end(); eter++)
	{
		if (primal.is_boundary(*eter))
		{
			PMyMesh::Point p(0, 0, 0);
			primal.property(edgev, *eter) = dual->add_vertex(p);
		}
	}
	/*for (PMyMesh::EdgeIter eter=primal.edges_begin();eter!=primal.edges_end();eter++)
	{

	}*/
	vector<VertexHandle> faces_handles;
	for (TMyMesh::VertexIter vter = primal.vertices_begin(); vter != primal.vertices_end(); vter++)
	{
		if ((*vter).idx()==(primal.n_vertices()-1))
		{
			continue;
		}
		if (!primal.is_boundary(*vter))
		{
			faces_handles.clear();
			for (TMyMesh::VertexFaceIter vfter = primal.vf_iter(*vter); vfter.is_valid(); vfter++)
			{
				faces_handles.push_back(primal.property(primalToDual, *vfter));
			}
			dual->add_face(faces_handles);
		}
		else
		{
			faces_handles.clear();
			for (TMyMesh::VertexFaceIter vfter = primal.vf_iter(*vter); vfter.is_valid(); vfter++)
			{
				faces_handles.push_back(primal.property(primalToDual, *vfter));
			}
			TMyMesh::VertexOHalfedgeCCWIter hter = primal.voh_ccwbegin(*vter);
			hter++;
			if (!primal.is_boundary(primal.edge_handle(*hter)))
			{
				printf_s("meshdual cuoiwu\r\n");
				int mytest;
				scanf_s("%d", &mytest);
			}
			faces_handles.push_back(primal.property(edgev, primal.edge_handle(*hter)));
			hter--;
			if (!primal.is_boundary(primal.edge_handle(*hter)))
			{
				printf_s("meshdual cuoiwu\r\n");
				int mytest;
				scanf_s("%d", &mytest);
			}
			faces_handles.push_back(primal.property(edgev, primal.edge_handle(*hter)));
			dual->add_face(faces_handles);
		}
	}
	primal.remove_property(primalToDual);
	return dual;
}

class Omt
{
public:
	Omt(TMyMesh &me, igl::viewer::Viewer &vie) :mesh(me), viewer(vie)
	{

		printf_s("create omt\r\n");
	}
	~Omt()
	{}
	void set_mesh1()
	{
		mesh1 = set_meshdual(mesh);
		
	}
	
	Eigen::VectorXd gradient_of_face_new(FaceHandle f)
	{
		double f_star[3];
		TMyMesh::Point p[3];
		Eigen::VectorXd temp_v;
		Eigen::Matrix2d temp_m;
		HalfedgeHandle he = mesh.halfedge_handle(f);
		if (mesh.from_vertex_handle(he).idx()==(mesh.n_vertices()-1)||mesh.to_vertex_handle(he).idx()==(mesh.n_vertices()-1)||mesh.to_vertex_handle(mesh.next_halfedge_handle(he)).idx()==(mesh.n_vertices()-1))
		{
			if (mesh.from_vertex_handle(he).idx() == (mesh.n_vertices() - 1))
			{
				he = mesh.next_halfedge_handle(he);
			}
			else if (mesh.to_vertex_handle(he).idx() == (mesh.n_vertices() - 1))
			{
				he = mesh.prev_halfedge_handle(he);
			}
			else
			{
			}
			
			Eigen::Vector3d temp_v1, temp_v2, temp_v3;
			//Eigen::VectorXd temp_v;
			temp_v1.setZero();
			temp_v2.setZero();
			temp_v3.setZero();
			temp_v = gradient_of_face_new(mesh.face_handle(mesh.opposite_halfedge_handle(he)));
			temp_v1.coeffRef(0) = temp_v.coeff(0);
			temp_v1.coeffRef(1) = temp_v.coeff(1);
			VertexHandle v1, v2;
			v1 = mesh.to_vertex_handle(he);
			v2 = mesh.from_vertex_handle(he);
			temp_v2.coeffRef(0) = mesh.point(v1)[0] - mesh.point(v2)[0];
			temp_v2.coeffRef(1) = mesh.point(v1)[1] - mesh.point(v2)[1];

			//temp_v2.normalize();
			temp_v3 = (temp_v1.transpose()*temp_v2.normalized()).coeff(0)*temp_v2.normalized();
			return (temp_v3 + (temp_v1 - temp_v3).normalized()*sqrt(1 - temp_v3.squaredNorm()));
			
		}
		else {
			temp_v.resize(2,1);
			int i = 0;
			for (TMyMesh::FaceVertexIter fvter = mesh.fv_begin(f); fvter.is_valid(); fvter++)
			{
				p[i] = mesh.point(*fvter);
				//f_star[i] = mesh.data(*fvter).f_star;
				f_star[i] = F_STAR.coeff((*fvter).idx());
				i++;
			}
			for (i = 0; i<2; i++)
			{
				for (int j = 0; j<2; j++)
				{
					temp_m.coeffRef(i, j) = p[i + 1][j] - p[i][j];
				}
				temp_v.coeffRef(i) = f_star[i + 1] - f_star[i];
			}
		}
		
		return temp_m.inverse()*temp_v;
	
	}
	
	double edge_delaunay_new(EdgeHandle e)
	{
		HalfedgeHandle he = mesh.s_halfedge_handle(e, 0);
		Eigen::Vector3d temp_v1, temp_v2, temp_v3;
		Eigen::VectorXd temp_v;
		if (mesh.is_boundary(e))
		{
			printf_s("edge_delaunay\r\n");
		}
		else if (mesh.from_vertex_handle(he).idx() != (mesh.n_vertices() - 1) && mesh.to_vertex_handle(he).idx() != (mesh.n_vertices() - 1))
		{

			temp_v3.coeffRef(0) = mesh.point(mesh.to_vertex_handle(he))[0] - mesh.point(mesh.from_vertex_handle(he))[0];
			temp_v3.coeffRef(1) = mesh.point(mesh.to_vertex_handle(he))[1] - mesh.point(mesh.from_vertex_handle(he))[1];
			temp_v3.coeffRef(2) = 0;
			temp_v = gradient_of_face_new(mesh.face_handle(he));
			temp_v1.coeffRef(0) = temp_v.coeff(0);
			temp_v1.coeffRef(1) = temp_v.coeff(1);
			temp_v1.coeffRef(2) = 0;
			temp_v = gradient_of_face_new(mesh.face_handle(mesh.opposite_halfedge_handle(he)));
			temp_v2.coeffRef(0) = temp_v.coeff(0);
			temp_v2.coeffRef(1) = temp_v.coeff(1);
			temp_v2.coeffRef(2) = 0;
			return temp_v3.cross(temp_v1 - temp_v2).coeff(2);
		}
		else
		{
			if (mesh.from_vertex_handle(he).idx() == (mesh.n_vertices() - 1))
			{
				he = mesh.opposite_halfedge_handle(he);
			}
			else if (mesh.to_vertex_handle(he).idx() == (mesh.n_vertices() - 1))
			{

			}
			else
			{
				printf_s("edge_delaunay_new\r\n");
				return 0;
			}
			temp_v = gradient_of_face_new(mesh.face_handle(he));
			temp_v1.coeffRef(0) = temp_v.coeff(0);
			temp_v1.coeffRef(1) = temp_v.coeff(1);
			temp_v1.coeffRef(2) = 0;
			temp_v = gradient_of_face_new(mesh.face_handle(mesh.opposite_halfedge_handle(he)));
			temp_v2.coeffRef(0) = temp_v.coeff(0);
			temp_v2.coeffRef(1) = temp_v.coeff(1);
			temp_v2.coeffRef(2) = 0;
			temp_v3.coeffRef(0) = temp_v1.coeff(0) + temp_v2.coeff(0);
			temp_v3.coeffRef(1) = temp_v1.coeff(1) + temp_v2.coeff(1);
			temp_v3.coeffRef(2) = 0;
			return temp_v3.cross(temp_v1 - temp_v2).coeff(2);
		}
	}
	
	int  adjust_new()
	{
		//printf_s("mesh.n_edges : %d\r\n", mesh.n_edges());
		int  temp_b = 0;
		for (TMyMesh::VertexIter vter = mesh.vertices_begin(); vter != mesh.vertices_end(); vter++)
		{
			for (TMyMesh::VertexOHalfedgeCCWIter vhter = mesh.voh_ccwbegin(*vter); vhter.is_valid(); vhter++)
			{
				HalfedgeHandle he = mesh.next_halfedge_handle(*vhter);
				if (mesh.is_boundary(mesh.s_edge_handle(he)))
				{

				}
				else if (edge_delaunay_new(mesh.s_edge_handle(he))<0)
				{
					//printf_s("edge_delaunay %lf\r\n", edge_delaunay(mesh.s_edge_handle(he)));

					if (mesh.is_flip_ok(mesh.s_edge_handle(he)))
					{
						//printf_s("adjust test ok\r\n");
						temp_b = 1;
						mesh.flip(mesh.s_edge_handle(he));
					}
					else
					{
					
						//printf_s("adjust test\r\n");
						temp_b = -1;
						return temp_b;
					}
				}
				else
				{
				}
				//FaceHandle f;
				//Eigen::Vector3d temp_v1, temp_v2;
				//f = mesh.face_handle(*vhter);
				//if (mesh.is_valid_handle(f))
				//{
				//	if (!mesh.is_boundary(mesh.s_edge_handle(mesh.next_halfedge_handle(*vhter))))
				//    {
				//		temp_v1 = mesh.data(f).v;
				//		f = mesh.face_handle(mesh.opposite_halfedge_handle(mesh.next_halfedge_handle(*vhter)));
				//		temp_v2 = mesh.data(f).v;
				//		Eigen::Vector3d temp_v3; TMyMesh::HalfedgeHandle he = mesh.next_halfedge_handle(*vhter);
				//		temp_v3.coeffRef(0) = mesh.point(mesh.to_vertex_handle(he))[0]-mesh.point(mesh.from_vertex_handle(he))[0];
				//		temp_v3.coeffRef(1) = mesh.point(mesh.to_vertex_handle(he))[1] - mesh.point(mesh.from_vertex_handle(he))[1];
				//		temp_v3.coeffRef(2) =0;

				//		if (temp_v3.cross(temp_v1 - temp_v2).coeff(2)<1e-16)
				//		{
				//			printf_s("adjust test %f\r\n", temp_v3.cross(temp_v1 - temp_v2).coeff(2));
				//			//TMyMesh test;
				//			if (mesh.is_flip_ok(mesh.s_edge_handle(he)))
				//			{
				//				mesh.flip(mesh.s_edge_handle(he));
				//			}
				//			else
				//			{
				//				printf_s("adjust test\r\n");
				//			}
				//		}
				//	}
				//}
				//else
				//{
				//}

			}
		}
		return temp_b;
		printf_s("%d\r\n", mesh.n_edges());
	}
	
	void update_dual_point_new()
	{
		
		int temp_i = 0;
		PMyMesh::VertexIter vter = mesh1->vertices_begin();
		for (TMyMesh::FaceIter fter = mesh.faces_begin(); fter != mesh.faces_end(); fter++)
		{
			Eigen::VectorXd v1 = gradient_of_face_new(*fter);
			/*mesh.data(*fter).v.coeffRef(0) = v1.coeff(0);
			mesh.data(*fter).v.coeffRef(1) = v1.coeff(1);
			mesh.data(*fter).v.coeffRef(2) = 0;*/
			mesh1->point(*vter)[0]= v1.coeff(0);
			mesh1->point(*vter)[1] = v1.coeff(1);
			mesh1->point(*vter)[2] = 0;
			dual_V.coeffRef(temp_i, 0) = v1.coeff(0);
			dual_V.coeffRef(temp_i, 1) = v1.coeff(1);
			dual_V.coeffRef(temp_i, 2) = 0;
			temp_i++;
			vter++;
		}
	
	}
	
	
	void init()
	{
		set_mesh1();
		A.resize(mesh.n_vertices(), 1);
		W.resize(mesh.n_vertices(), 1);
		F_STAR.resize(mesh.n_vertices(), 1);
		Fi.resize(mesh.n_vertices(), 1);
		H.resize(mesh.n_vertices()-1,mesh.n_vertices()-1);
		dual_V.resize(mesh1->n_vertices(), 3);
	
		int i = 0;
		for (TMyMesh::VertexIter vter = mesh.vertices_begin(); vter != mesh.vertices_end(); vter++)
		{
				F_STAR.coeffRef(i) = mesh.point(*vter).sqrnorm() / 2.0;
				//mesh.data(*vter).f_star = F_STAR.coeffRef(i);
				i++;
		}
		//增加额外点，其点代表圆柱面
		TMyMesh::Point p(0, 0, -1);
		VertexHandle v=mesh.add_vertex(p);
		printf_s("vid %d,n_vertices-1 %d\r\n", v.idx(),mesh.n_vertices()-1);
		vector<VertexHandle> face_handles;
		list<vector<VertexHandle>> lface_handles;
		
		for (TMyMesh::EdgeIter eter=mesh.edges_begin();eter!=mesh.edges_end();eter++)
		{
			if (mesh.is_boundary(*eter))
			{
				face_handles.clear();
				HalfedgeHandle he;
				if (mesh.is_valid_handle(mesh.face_handle(mesh.s_halfedge_handle(*eter,0))))
				{
					he = mesh.s_halfedge_handle(*eter,1);
				}
				else
				{
					he = mesh.s_halfedge_handle(*eter,0);
				}
				face_handles.push_back(mesh.from_vertex_handle(he));
				face_handles.push_back(mesh.to_vertex_handle(he));
				face_handles.push_back(v);
				
				lface_handles.push_back(face_handles);
			}
			else
			{
			}
		}
		printf_s("lface_handle.seize %d\r\n",lface_handles.size());
		for (list<vector<VertexHandle>>::iterator lfter=lface_handles.begin();lfter!=lface_handles.end();lfter++)
		{
			mesh.add_face(*lfter);
		}
		
	}
	Eigen::VectorXd oriarea(vector<Eigen::VectorXd> vv)
	{
		int n = vv.size();
		if (n<3)
		{
			printf_s("oriarea cuowu\r\n");
			int mytest;
			scanf_s("%d", &mytest);
		}
		Eigen::Vector3d temp = vv[0], sum;
		sum.setZero();
		for (int i = 1; i<n; i++)
		{
			Eigen::Vector3d temp1, temp2;
			temp1 = vv[i] - temp;
			temp2 = vv[(i + 1)%n] - vv[i];
			sum += temp1.cross(temp2);
		}
		return sum / 2.0;
		//sum.cros

	}
	void area_dual()
	{
		//printf_s("mesh1face %d\r\n",mesh1->n_faces());
		/*for (TMyMesh::VertexIter vter=mesh.vertices_begin();vter!=mesh.vertices_end();vter++)
		{i
		}*/
		int i = 0;
		vector<Eigen::VectorXd>vv;
		for (PMyMesh::FaceIter fter = mesh1->faces_begin(); fter != mesh1->faces_end(); fter++)
		{
			vv.clear();
			for (PMyMesh::FaceVertexIter fvter=mesh1->fv_begin(*fter);fvter.is_valid();fvter++)
			{
				Eigen::Vector3d temp_v;
				temp_v.coeffRef(0) = dual_V.coeff((*fvter).idx(), 0);
				temp_v.coeffRef(1) = dual_V.coeff((*fvter).idx(), 1);
				temp_v.coeffRef(2) = 0;
				vv.push_back(temp_v);
			}
			A.coeffRef(i,0)=oriarea(vv).norm();
			
			//printf_s("A(i) %lf\r\n", A.coeffRef(i,0));(Eigen::MatrixXd::Ones(1, A.rows())*A).coeff(0)
			i++;
		}
		A *= (acos(-1)/A.sum());
		printf_s("A sum %lf\r\n", A.sum() );

	}
	
	void compute_hessian()
	{
		printf_s("hessian\r\n");
		H.setZero();
		for (TMyMesh::EdgeIter eter=mesh.edges_begin();eter!=mesh.edges_end();eter++)
		{
			HalfedgeHandle he = mesh.s_halfedge_handle(*eter,0);
			VertexHandle v1= mesh.to_vertex_handle(he), v2= mesh.from_vertex_handle(he);
			if (v1.idx()==(mesh.n_vertices()-1)|| v2.idx() == (mesh.n_vertices() - 1))
			{
				continue;
			}
			else
			{
				FaceHandle f = mesh.face_handle(mesh.opposite_halfedge_handle(he));
				double tempd = -(dual_V.row(mesh.face_handle(he).idx()) - dual_V.row(f.idx())).norm();
				tempd /= mesh.calc_edge_length(*eter);
				/*H.coeffRef(v1.idx(), v2.idx()) = tempd;
				H.coeffRef(v2.idx(), v1.idx()) = tempd;
				H.coeffRef(v1.idx(), v1.idx()) -= tempd;
				H.coeffRef(v2.idx(), v2.idx()) -= tempd;*/
				if (v1.idx()==(mesh.n_vertices()-2))
				{
					
					H.coeffRef(v2.idx(), v2.idx()) -= tempd;
				}
				else if(v2.idx()== (mesh.n_vertices() - 2))
				{
					H.coeffRef(v1.idx(), v1.idx()) -= tempd;
					
				}
				else
				{
					H.coeffRef(v1.idx(), v2.idx()) = tempd;
					H.coeffRef(v2.idx(), v1.idx()) = tempd;
					H.coeffRef(v1.idx(), v1.idx()) -= tempd;
					H.coeffRef(v2.idx(), v2.idx()) -= tempd;
				}
			}
		}
		printf_s("compressed begin\r\n");
		H.makeCompressed();
	
		printf_s("hessian end\r\n");
	}
	
	void iteration1()
	{
		double t = 0.25, e=1;
		
	   while (e>0.00000076)
	   {
		   if (e<0.0000009)
		   {
			   t = 0.255;
		   }
		   else if (e<0.0000008)
		   {
			   t = 0.27;
		   }
		   else if(e<0.00000079)
		   {
			   t = 0.28;
		   }


			//adjust_new();
			update_dual_point_new();
			area_dual();

			Eigen::VectorXd temp_v = A - W;

			F_STAR += t * temp_v;

			/*for (int i = 0; i<temp_v.rows(); i++)
			{
				mesh.data(mesh.vertex_handle(i)).f_star = F_STAR.coeff(i);
			}*/
			e = temp_v.norm();
			printf_s("energy: %.12lf\r\n", e);
		}
	  
	   OpenMesh::IO::Options opt;
	   OpenMesh::IO::write_mesh(*mesh1,"opt.off",opt,10);
	}
	void iteration2()
	{
		//A = A * (acos(-1)/A.norm());
		Eigen::VectorXd temp_v = (A - W), y;
		y.resize(F_STAR.rows(), 1);
		y.setZero();
		int temp_d;
		double t = 0.3, e = 1;
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>llt;
		//for(int j=0;j<14;j++)
		while (e>4e-11)
		{
			temp_d = adjust_new();
			//while(temp_d==-1) 
			//{
			//	//printf_s("temp_d%d\r\n",temp_d);
			//	F_STAR -= 0.005*t * y;
			//	/*for (int i = 0; i<A.rows(); i++)
			//	{
			//		mesh.data(mesh.vertex_handle(i)).f_star = F_STAR.coeff(i);
			//	}*/
			//	temp_d = adjust_new();
			//}
			
			if (temp_d != 0)
			{
				set_mesh1();
			}
			
			//adjust_new();
			update_dual_point_new();
			area_dual();
			compute_hessian();
			//A.head(A.rows()-1);
			
			temp_v = A - W;
			printf_s("temp_v sum %.16lf\n", temp_v.sum());
			e = temp_v.norm();
			printf_s("energy: %.15lf\r\n", e);
			
			llt.compute(H);
			y.head(F_STAR.rows()-1) = llt.solve(temp_v.head(F_STAR.rows()-1));
			F_STAR += t * y;
			
			printf_s("y.coeff %.10lf    %.10lf  average %.10lf\r\n", y.minCoeff(), (-1 * y).minCoeff(),y.sum()/double(y.rows()-1));
			printf_s("test %.15lf,,y.norm %.15lf\r\n",(-H*y.head(F_STAR.rows() - 1) + temp_v.head(F_STAR.rows() - 1)).norm(),y.norm());
			
		}
		OpenMesh::IO::Options opt;
		OpenMesh::IO::write_mesh(*mesh1, "opt.off", opt, 10);
		//printf_s("y.coeff %.10lf    %.10lf\r\n", y.minCoeff(),(-1*y).minCoeff());
			//mesh.data(mesh.vertex_handle(i)).f_star = F_STAR.coeff(i);
		

	}
	
	void drawaxis()
	{
		Eigen::RowVector3d v1, v2;
		v1.setZero(); v2.setZero();
		v2.coeffRef(2) = 1;
		viewer.data.add_edges(v1, v2, Eigen::RowVector3d(0.5, 1, 0));
		v2.coeffRef(2) = 0; v2.coeffRef(1) = 1;
		viewer.data.add_edges(v1, v2, Eigen::RowVector3d(0, 1, 0));
		for (int i = 0; i<100; i++)
		{
			v1 = Eigen::RowVector3d(cos(acos(-1)*i / 50.0), sin(acos(-1)*i / 50.0), 0);
			v2 = Eigen::RowVector3d(cos(acos(-1)*(i + 1) / 50.0), sin(acos(-1)*(i + 1) / 50.0), 0);
			viewer.data.add_edges(v1, v2, Eigen::RowVector3d(0.5, 1, 0));
		}
	}
	
	void drawedge()
	{
		printf_s("begin drawedge\r\n");
		Eigen::RowVector3d v1, v2; VertexHandle vv1, vv2;
		HalfedgeHandle he;
		for (PMyMesh::EdgeIter eter = mesh1->edges_begin(); eter != mesh1->edges_end(); eter++)
		{
			he = mesh1->s_halfedge_handle(*eter, 0);
			vv1 = mesh1->to_vertex_handle(he);
			vv2 = mesh1->from_vertex_handle(he);
			v1 = dual_V.row(vv1.idx());
			v2 = dual_V.row(vv2.idx());
			viewer.data.add_edges(v1, v2, Eigen::RowVector3d(0.5, 1, 0));
		}

	}
	
	
	Eigen::MatrixXd V, dual_V;
	Eigen::SparseMatrix<double> H;

	Eigen::MatrixXi F;
	Eigen::VectorXd A, W, F_STAR, Fi;
	PMyMesh *mesh1;
	double Tt = 0.0002, energy = 0;
private:
	igl::viewer::Viewer &viewer;
	TMyMesh & mesh;

};

