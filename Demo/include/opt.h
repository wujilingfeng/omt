
#ifndef LIBCELL_OPT
#define LIBCELL_OPT

#include<Mesh_IO/Mesh_IO.h>

#include "Arcroll.h"
#include<Algorithm/Algorithm.h>
#include<Sparse>
#include<Cholesky>
#include<SparseQR>
typedef struct Libcell_Opt
{
     Node* node_v;
    int v_size,rows;
    template_v*  v0;
    Mesh* mesh;
    double sum_meas;
    double**H;
    double* Wi,*Ai;

}L_Opt;
void Libcell_Opt_init(L_Opt* op);
//在每点prop给个数值
void libcell_compute_dual_point(template_m*,template_c* c_);

bool convex_adjust(template_m*,template_f*);

void iteration(L_Opt*myopt);
//返回边界不是凸的顶点
#endif
