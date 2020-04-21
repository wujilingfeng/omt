#include<Mesh_view/Mesh_viewer_interpreter.h>
#include "include/opt.h"
#include<tool/libcell_tools_view.h>
#include<iostream>
typedef struct myanimation
{
    Mesh* mesh;
    Mesh_viewer_faces*mf;
    Mesh_viewer_opengl_Interpreter* moi;
    int i;
    double **v;
}myani;
//检测局部凸性，检测局部形成的闭凸包的边界的方向与原方向是否相反
void test3D();
void test_camera_and_intera(Mesh_viewer_world*);
void test_delauny();
void test_convex();
void once(Mesh_viewer_Intera* );

void convex_init(Mesh*,double **,int);
void test_area();
void test_opt();
int main()
{
   /* Mesh mesh;
    A_Mesh_init(&mesh);
    _ReadCell_(&mesh,"delauny_circle.cell");
    Node*v=mesh.isolate_vertices(&mesh);
    printf("size:%d\n",node_size(v));*/
    /*v=(Node*)(v->Next); v=(Node*)(v->Next);
    v=(Node*)(v->Next);

    if(node_size(v)!=0)
    {
        printf("buhao:%d\n",((template_v*)(v->value))->id);
    }
    for(auto iter=mesh.cells.begin();iter!=mesh.cells.end();iter++)
    {
        for(auto it=mesh.cv_begin(&mesh,*(iter->second));it!=mesh.cv_end(&mesh,*(iter->second));it++)
        {
            if((*it).id==156)
            {   
                printf("zhenbuah\n");
            }
        }
    }*/
//    Eigen::MatrixXd A(24,4);
/*    A<<0.070520521200000005346275600004,0.092321399999999997909583271394,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.113275185400000005797593871648,0.064418742000000001057813392435,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.082658631999999995776384764667,0.106561249900000001855637776771,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.080065807200000005394713298301,0.106625704900000006780658168282,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.082562390099999993786639151949,0.089196737999999997636990656247,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.114646251399999996034573257475,0.054330086399999998059140438045,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.087457438100000006730994073223,-0.016363862800000000125022481257,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.073725474499999998911903276166,-0.006239110800000000241205544427,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.068192995100000000641138342417,0.016326049100000000063870331246,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.086984531200000006490036241757,0.043811974500000003263089354277,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.049885169399999998440797810417,0.021983758900000000685581724724,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.032954001099999997992373579336,0.042383597000000002019692857402,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.026021966099999998406344658974,0.056822831800000001001293981062,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.016340481000000000449956516491,0.053546509600000000150821222178,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.019678191700000000063930372107,0.046811146400000003264807446612,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.030237093499999999346439238934,0.009181909199999999451891063984,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.029913883499999998488716457246,0.007038980700000000041705039422,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.020388770800000000776863728902,0.004971218000000000337978978138,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.014707106799999999577344134138,0.023433714700000000258350851823,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.013670972499999999721920218576,0.054224362200000000844291037083,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.038907190799999998020997082904,0.020418590199999999551527807284,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.042046912300000002793431974624,-0.005748506599999999785999271751,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.041157782400000002265283427505,-0.004181704900000000370285935247,-0.214748364800000007823754799574,0.000000000000000000000000000000,
0.036759981999999996604699958880,-0.003287679200000000016845813633,-0.214748364800000007823754799574,0.000000000000000000000000000000;*/
  /*      A<<0.0705205212,0.0923213999,-0.2147483648,0.0000000000,
0.1132751854,0.0644187420,-0.2147483648,0.0000000000,
0.0826586319,0.1065612499,-0.2147483648,0.0000000000,
0.0800658072,0.1066257049,-0.2147483648,0.0000000000,
0.0825623900,0.0891967379,-0.2147483648,0.0000000000,
0.1146462513,0.0543300863,-0.2147483648,0.0000000000,
0.0874574381,-0.0163638628,-0.2147483648,0.0000000000,
0.0737254744,-0.0062391108,-0.2147483648,0.0000000000,
0.0681929951,0.0163260491,-0.2147483648,0.0000000000,
0.0869845312,0.0438119745,-0.2147483648,0.0000000000,
0.0498851693,0.0219837589,-0.2147483648,0.0000000000,
0.0329540010,0.0423835970,-0.2147483648,0.0000000000,
0.0260219660,0.0568228318,-0.2147483648,0.0000000000,
0.0163404810,0.0535465096,-0.2147483648,0.0000000000,
0.0196781917,0.0468111464,-0.2147483648,0.0000000000,
0.0302370934,0.0091819091,-0.2147483648,0.0000000000,
0.0299138834,0.0070389807,-0.2147483648,0.0000000000,
0.0203887708,0.0049712180,-0.2147483648,0.0000000000,
0.0147071067,0.0234337147,-0.2147483648,0.0000000000,
0.0136709724,0.0542243622,-0.2147483648,0.0000000000,
0.0389071907,0.0204185901,-0.2147483648,0.0000000000,
0.0420469123,-0.0057485065,-0.2147483648,0.0000000000,
0.0411577824,-0.0041817049,-0.2147483648,0.0000000000,
0.0367599819,-0.0032876792,-0.2147483648,0.0000000000;
    //std::cout<<A<<std::endl;
    for(int i=0;i<3;i++)
    {
      //  printf("%lf ",A.data()[i]);
    }
    Mesh mesh;
    A_Mesh_init(&mesh);
    double **VV=(double**)malloc(sizeof(double*)*24);
    for(int i=0;i<24;i++)
    {
        VV[i]=(double*)malloc(sizeof(double)*4);
        for(int j=0;j<4;j++)
        {
            VV[i][j]=A.data()[j*24+i];
            //int a=VV[i][j]*10000000000;
            //VV[i][j]=a/10000000000.0;
            printf("%lf ",VV[i][j]);
        }
        printf("\n");
    }
    compute_convex_area(VV,24,4,3);
*/
  // double t=0.2345;
  //  int a=t*100;
   // printf("%d  \n",a);
    test_opt();
    //test_area();
//    test3D(); 
    //test_convex();
//    test_delauny();
/*    Mesh mesh;
    Mesh_init(&mesh);
    double **v=(double**)malloc(sizeof(double*)*8);
    for(int i=0;i<8;i++)
    {
        v[i]=(double*)malloc(sizeof(double)*3);

    }
    v[0][0]=-0.5;v[0][1]=-0.5;v[0][2]=0.5;
    v[1][0]=0.5;v[1][1]=-0.5;v[1][2]=0.5;
    v[2][0]=-0.5;v[2][1]=0.5;v[2][2]=0.5;
    v[3][0]=-0.5;v[3][1]=0.5;v[3][2]=-0.5;
    v[4][0]=0.5;v[4][1]=0.5;v[4][2]=0.5;
    v[5][0]=0.5;v[5][1]=0.5;v[5][2]=-0.5;
    v[6][0]=-0.5;v[6][1]=-0.5;v[6][2]=-0.5;
    v[7][0]=0.5;v[7][1]=-0.5;v[7][2]=-0.5;
    from_v_createdelauny_simplex(&mesh,v,8,3);
    for(int i=0;i<8;i++)
    {
        free(v[i]);
    }
    free(v);
*/
   /* Mesh mesh;
    Mesh_init(&mesh);
    _ReadCell_(&mesh,"bunny.cell");
    auto iter=mesh.cells.begin();
    for(auto iter1=mesh.chf_begin(&mesh,*(iter->second));iter1!=mesh.chf_end(&mesh,*(iter->second));iter1++)
    {
        for(auto iter2=mesh.hfv_begin(&mesh,*(iter1));iter2!=mesh.hfv_end(&mesh,*(iter1));iter2++)
        {
            printf("%d  ",(*iter2).id); 
        }
        printf("\n");
    
    }*/
    //moi.render(&moi);
    //float*v
    //printf("%d\n",mesh.num_v(&mesh));    
       return 0;
}
static void init_opt1(L_Opt*myopt)
{
    Mesh* mesh=myopt->mesh;
    myopt->v_size=mesh->num_v(mesh);
    myopt->Ai=(double*)malloc(sizeof(double)*myopt->v_size);
    for(int i=0;i<myopt->v_size;i++)
    {
        myopt->Ai[i]=1.0/(double)(myopt->v_size);
    }
    printf("num_v:%d\n",mesh->num_v(mesh));
    double temp_sum=0;
    //double *v=NULL;
    int *id=NULL;
    for(auto iter=mesh->vertices.begin();iter!=mesh->vertices.end();iter++)
    {
        temp_sum=0;
        double* temp_v=(double*)malloc(sizeof(double)*(iter->second->point_size+1));
        for(int j=0;j<iter->second->point_size;j++)
        {
            temp_v[j]=iter->second->point[j];
            temp_sum+=iter->second->point[j]*iter->second->point[j];
        }
        temp_v[iter->second->point_size]=temp_sum/2.0;
        iter->second->point_size++;
        free(iter->second->point);
        iter->second->point=temp_v;
        //v=(double*)malloc(sizeof(double));
        //iter->second->prop=(void*)v;
        id=(int*)malloc(sizeof(int));
        iter->second->prop1=(void*)id;
    }
    for(auto iter=mesh->cells.begin();iter!=mesh->cells.end();)
    {
        mesh->delete_cell(mesh,*((iter++)->second),true); 
    }
    template_v*v0=NULL;
    double *temp_v=(double*)malloc(sizeof(double)*(mesh->dimension+1));
    memset(temp_v,0,sizeof(double)*(mesh->dimension));
    temp_v[mesh->dimension]=10000.0;
    v0=mesh->create_vertexv(mesh,temp_v,mesh->dimension+1);
    myopt->v0=v0;
}
static void init_opt(Mesh*mesh)
{
    double temp_sum;
    double *temp_v=NULL;
    for(auto iter=mesh->vertices.begin();iter!=mesh->vertices.end();iter++)
    {
        temp_sum=0;
        double* f_i=(double*)malloc(sizeof(double));
        iter->second->prop=(void*)f_i;
        temp_v=(double*)malloc(sizeof(double)*(iter->second->point_size+1));
        for(int j=0;j<iter->second->point_size;j++)
        {
          temp_v[j]=iter->second->point[j];
            temp_sum+=iter->second->point[j]*iter->second->point[j];
        }

        *f_i=temp_sum/2.0;
        if(iter->second->id==500)
        {
            *f_i=0.0;
        }
        temp_v[iter->second->point_size]=*f_i;
        iter->second->point_size++;
        free(iter->second->point);
        iter->second->point=temp_v;
    }
    template_v*v0=NULL;
    temp_v=(double*)malloc(sizeof(double)*(mesh->dimension+1));
    memset(temp_v,0,sizeof(double)*(mesh->dimension));
    temp_v[mesh->dimension]=10000.0;
    v0=mesh->create_vertexv(mesh,temp_v,mesh->dimension+1);
    free(temp_v);
    template_v** temp_v1=(template_v**)malloc(sizeof(template_v*)*(mesh->dimension+1));
    mesh->external_cell_init_(mesh);
    Node* n_it=mesh->external_cell.halffaces;
    template_hf* hf0=NULL;
    while(n_it!=NULL)
    {
        hf0=(template_hf*)(n_it->value);
        for(int i=0;i<hf0->vertices_size;i++)
        {
            temp_v1[i]=(template_v*)(hf0->vertices[i]);
        }
        temp_v1[mesh->dimension]=v0;
        mesh->create_cellv(mesh,temp_v1,mesh->dimension+1);
        n_it=(Node*)(n_it->Next);
    } 
}
static void cut_out(Mesh* mesh)
{   
    double temp_sum=0;
    double*p=NULL;
    Node* node_c=NULL;
    for(auto iter=mesh->cells.begin();iter!=mesh->cells.end();iter++)
    {
        temp_sum=0;

        p=(double*)(iter->second->prop);
        for(int j=0;j<mesh->dimension;j++)
        {
            temp_sum+=p[j]*p[j];
        }
        if(temp_sum>1)
        {
            node_c=node_overlying(node_c,iter->second);
        //    mesh->delete_cell(mesh,*((iter++)->second),true); 
        } 
    }
    while(node_c!=NULL)
    {
    
        mesh->delete_cell(mesh,*((template_c*)(node_c->value)),true);
        node_c=(Node*)(node_c->Next);
    }

}
static void improve_topology_from_mesh(L_Opt*myopt)
{
    Mesh* mesh=myopt->mesh;
    template_v* v0=myopt->v0;
    mesh->external_cell_init_(mesh);
    Node* n_it=mesh->external_cell.halffaces;
    int dim=mesh->dimension;

    double * temp_v=NULL;
    template_hf* hf0=NULL,*hf1=NULL;
    Eigen::MatrixXd A,temp_re,temp_b;
    template_v *v1=NULL,*vi=NULL;
    double a,b,c;
    double *p,*p1;
    A.resize(dim-1,dim);
    while(n_it!=NULL)
    {

        temp_v=(double*)malloc(sizeof(double)*dim);
        hf0=(template_hf*)(n_it->value);
        hf0->prop=(void*)temp_v;
        hf1=mesh->s_opposite_halfface(hf0);
        v1=(template_v*)(hf0->vertices[0]);
        for(int i=1;i<hf0->vertices_size;i++)
        {
            vi=(template_v*)(hf0->vertices[i]);
            for(int j=0;j<vi->point_size-1;j++)
            {
            
                A.coeffRef(i-1,j)=vi->point[j]-v1->point[j];
            }
        }

        p=temp_normal_vector(A);p1=(double*)(((template_c*)(hf1->cell))->prop);
        a=0;b=0;c=0;
        for(int i=0;i<dim;i++)
        {
           a+=p[i]*p[i];
           b+=(2*p[i]*p1[i]);
           c+=p1[i]*p1[i];
        }

        c=c-1;
        double t=(-b+sqrt(b*b-4*a*c))/(2*a);
        for(int i=0;i<dim;i++)
        {
           temp_v[i]=p1[i]+p[i]*t;
        }
        n_it=(Node*)(n_it->Next);
    }
    n_it=mesh->external_cell.halffaces;
    template_c* c_=NULL;
    template_v**temp_v1=(template_v**)malloc(sizeof(template_v*)*(mesh->dimension+1));
    while(n_it!=NULL)
    {
        hf0=(template_hf*)(n_it->value);
        for(int j=0;j<hf0->vertices_size;j++)
        {
            temp_v1[j]=(template_v*)(hf0->vertices[j]);
        }
        temp_v1[mesh->dimension]=v0;
        c_=mesh->create_cellv(mesh,temp_v1,mesh->dimension+1);
        c_->prop=hf0->prop;
        
        n_it=(Node*)(n_it->Next);
    }
}
static void draw_points(Mesh_viewer_points* mp,Mesh* mesh)
{

    mp->Data_rows=mesh->num_c(mesh);
    mp->Data=(float*)malloc(sizeof(float)*3*mp->Data_rows);
    memset(mp->Data,0,sizeof(float)*3*mp->Data_rows);


}
static void draw_edges(Mesh_viewer_edges* me,Mesh* mesh)
{
    me->Data_index_rows=mesh->num_f(mesh);
    
    me->Data_rows=me->Data_index_rows*2;
    me->Data=(float*)malloc(sizeof(float)*3*me->Data_rows);
    me->Data_index=(unsigned int*)malloc(sizeof(unsigned int)*2*me->Data_index_rows);
    memset(me->Data_index,0,sizeof(unsigned int)*2*me->Data_index_rows);
    me->color_rows=me->Data_index_rows;
    float color[]={1.0,0.0,0.0};
    me->set_color(me,color);
    double *p;
    template_hf* hf0=NULL,*hf1=NULL;
    unsigned int i=0;
    for(auto iter=mesh->faces.begin();iter!=mesh->faces.end();iter++)
    {
        hf0=iter->second->halffaces[0];
        hf1=iter->second->halffaces[1];

            me->Data_index[i]=i;
        if(hf0->cell==NULL)
        {
            p=(double*)(hf0->prop);
        }
        else
        {
            p=(double*)(((template_c*)(hf0->cell))->prop);
        }
        for(int j=0;j<mesh->dimension;j++)
        {
            me->Data[i*3+j]=p[j];  
        }
        i++;
        me->Data_index[i]=i;

        if(hf1->cell==NULL)
        {
             p=(double*)(hf1->prop);
        }
        else
        {
            p=(double*)(((template_c*)(hf1->cell))->prop); 
        }
        for(int j=0;j<mesh->dimension;j++)
        {
            me->Data[i*3+j]=p[j];  
        }
        i++;
    }

}
static void test_derivate(L_Opt *myopt)
{
    int rows=node_size(myopt->node_v)-1;

    template_v** V=(template_v**)malloc(sizeof(template_v*)*rows);
    Node* n_it=myopt->node_v;
    for(int i=0;i<rows;i++)
    {
        V[i]=(template_v*)(n_it->value);
        n_it=(Node*)(n_it->Next);
    }
    double sum=0;
}
static void iteration_end(L_Opt*myopt)
{
    Mesh*mesh=myopt->mesh;
    for(auto iter=mesh->cells.begin();iter!=mesh->cells.end();)
    {
        free(iter->second->prop);
        mesh->delete_cell(mesh,*((iter++)->second),true); 
    }

    mesh_createconvex(mesh);
    printf("numc:%d\n",mesh->num_c(mesh));
    for(auto iter=mesh->vc_begin(mesh,*(myopt->v0));iter!=mesh->vc_end(mesh,*(myopt->v0));iter++)
    {
        mesh->delete_cell(mesh,*iter,true); 
    }
    for(auto iter=mesh->cells.begin();iter!=mesh->cells.end();iter++)
    {
        libcell_compute_dual_point(mesh,iter->second);
    }
    cut_out(mesh);
    improve_topology_from_mesh(myopt);

}
void test_opt()
{
    L_Opt myopt;
    Libcell_Opt_init(&myopt);
    Mesh mesh;
    A_Mesh_init(&mesh);
    myopt.mesh=&mesh;
    _ReadCell_(&mesh,"delauny_sphere.cell");
    //init_opt(&mesh);
    init_opt1(&myopt);
    printf("num_v:%d\n",mesh.num_v(&mesh));
    ///kaishi
    mesh_createconvex(&mesh);
    printf("numc:%d\n",mesh.num_c(&mesh));
    for(auto iter=mesh.vc_begin(&mesh,*(myopt.v0));iter!=mesh.vc_end(&mesh,*(myopt.v0));iter++)
    {
        mesh.delete_cell(&mesh,*iter,true); 
    }
    for(auto iter=mesh.cells.begin();iter!=mesh.cells.end();iter++)
    {
        libcell_compute_dual_point(&mesh,iter->second);
    }
    cut_out(&mesh);
    improve_topology_from_mesh(&myopt);
    printf("numc11: %d\n",mesh.num_v(&mesh));
    /*for(int i=0;i<330;i++)
    {
        iteration(&myopt);
        if(i==120)
        {
            break;
        }
        iteration_end(&myopt);
    }*/
    
      /* f0=(template_f*)(node_f->value);
    flag=convex_adjust(&mesh,f0);
    free_node(node_f);
    node_f=NULL;*/

   /* if(flag)
    {
        printf("chengong\n");
    }*/
   /* libcell_compute_dual_point(&mesh);*/
    Mesh_viewer_world mw;
    Mesh_viewer_world_init(&mw);
    char ch[]="faces";
    Node* n=mw.create_something(&mw,ch);
    Mesh_viewer_something*ms=(Mesh_viewer_something*)(n->value);
    Mesh_viewer_faces *mf=(Mesh_viewer_faces*)(ms->evolution);
    mf->color_rows=mesh.num_c(&mesh);

    float color[]={1.0,1.0,0.0};
    //mf->set_color(mf,color);
    mf->random_color(mf);
    get_data_from_2dim_cell(&mesh,&(mf->Data),&(mf->Data_index));
    mf->Data_rows=mesh.num_v(&mesh);
    mf->Data_index_rows=mesh.num_c(&mesh);
    mf->normal_rows=mf->Data_rows;
    mf->is_reversal_normal=1;
    ms->disappear=1;

    free_node(n);
    /*char edges[]="edges,edges";
    n=mw.create_something(&mw,edges);
    ms=(Mesh_viewer_something*)(n->value);
    auto me=(Mesh_viewer_edges*)(ms->evolution);

    //printf("here\n");
    draw_edges(me,&mesh);

    n=(Node*)(n->Prev);
    ms=(Mesh_viewer_something*)(n->value);
    ms->disappear=1;
    me=(Mesh_viewer_edges*)(ms->evolution);
    me->Data_rows=2;me->Data_index_rows=1;
    me->Data=(float*)malloc(sizeof(float)*3*me->Data_rows);
    me->Data_index=(unsigned int*)malloc(sizeof(unsigned int)*2*me->Data_index_rows);
    me->Data[0]=1;me->Data[1]=0.0;me->Data[2]=0.0;
    me->Data[3]=0;me->Data[4]=1.0;me->Data[5]=0.0;
    me->Data_index[0]=0;me->Data_index[1]=1;
    me->color_rows=me->Data_index_rows;
    me->set_color(me,color);*/
    char points[]="points";
    n=mw.create_something(&mw,points);
    ms=(Mesh_viewer_something*)(n->value);
    auto mp=(Mesh_viewer_points*)(ms->evolution);
    mp->Data_rows=mesh.num_c(&mesh);
    mp->Data=(float*)malloc(sizeof(float)*3*mp->Data_rows);
    memset(mp->Data,0,sizeof(float)*3*mp->Data_rows);
    int i=0,cols=mesh.vertices.begin()->second->point_size;

    printf("cols:%d\n",cols);
    for(auto iter=mesh.cells.begin();iter!=mesh.cells.end();iter++)
    {
        //temp_sum=0
        if(iter->second->prop==NULL)
        {
            continue;
        }
        for(int j=0;j<cols-1;j++)
        {
            mp->Data[i*3+j]=((double*)(iter->second->prop))[j];
            //temp_sum+=mp->Data[i*3+j]*mp->Data[i*3+j];
        }
    //    mp->Data[i*3+2]=temp_sum/2.0;
        i++;
    }
    //ms->disappear=1;
    free_node(n);
    char Intera[]="Intera";
    n=mw.create_something(&mw,Intera);
    ms=(Mesh_viewer_something*)(n->value);
    Mesh_viewer_Intera* mi=(Mesh_viewer_Intera*)(ms->evolution);
    myani anima;
    anima.mesh=&mesh;anima.mf=mf;
    //mi->key_callback=once;

    mi->representation=(void*)(&anima);
    free_node(n);

    test_camera_and_intera(&mw);
    Mesh_viewer_opengl_Interpreter moi;
    Mesh_viewer_opengl_Interpreter_init(&moi);
    moi.world=&mw;
    anima.moi=&moi;
    moi.routine_show(&moi);
    Mesh_free(&mesh);
    
//    libcell_compute_dual_point(&mesh);



}
void test_area()
{
    Mesh mesh;
    A_Mesh_init(&mesh);
    double ** v=(double**)malloc(sizeof(double*)*4);
    for(int i=0;i<4;i++)
    {
        v[i]=(double*)malloc(sizeof(double)*3); 
    }
    v[0][0]=1;v[0][1]=1;v[0][2]=3;
    v[1][0]=-1;v[1][1]=1;v[1][2]=1;
    v[2][0]=-1;v[2][1]=-1;v[2][2]=-3;
    v[3][0]=2;v[3][1]=-1;v[3][2]=0;
    double area=compute_convex_area(v,4,3,2); 
    printf("%lf\n",area);
}
void test_delauny()
{
    srand((unsigned)time(0));
    Mesh mesh;
    A_Mesh_init(&mesh);
    double **v=(double**)malloc(sizeof(double*)*2000);
    for(int i=0;i<2000;i++)
    {
        v[i]=(double*)malloc(sizeof(double)*3);
    }
    for(int i=0;i<800;i++)
    {
        double r=1,delta=(rand()%2000)/1000.0-1,theta=(rand()%1000)/1000.0;
        //theta=0.5;
        v[i][0]=r*sin(theta*M_PI)*cos(delta*M_PI);
        v[i][1]=r*sin(theta*M_PI)*sin(delta*M_PI);
        v[i][2]=r*cos(theta*M_PI); 
    }

    for(int i=800;i<2000;i++)
    {
        double r=(rand()%1000)/1000.0,delta=(rand()%2000)/1000.0-1,theta=(rand()%1000)/1000.0;
        v[i][0]=r*sin(theta*M_PI)*cos(delta*M_PI);
        v[i][1]=r*sin(theta*M_PI)*sin(delta*M_PI);
        v[i][2]=r*cos(theta*M_PI);
          
    }
    from_v_createdelauny_simplex(&mesh,v,2000,3);
    for(auto it=mesh.vertices.begin();it!=mesh.vertices.end();it++)
    {
        it->second->point_size--; 
    }

    Node*re=mesh.isolate_vertices(&mesh);
    Node* n_it=re;
    //printf("re size:%d\n",node_size(re));
    while(re!=NULL)
    {
        mesh.delete_vertex(&mesh,*((template_v*)(re->value)),true);
        re=(Node*)(re->Next);
    }
    free_node(n_it);
    _WriteCell_(&mesh,"delauny_sphere2.cell");
    


    //Mesh_viewer_world mw;
    //Mesh_viewer_world_init(&mw);
    /*char ch[]="faces";
    Node* n=mw.create_something(&mw,ch);
    Mesh_viewer_something*ms=(Mesh_viewer_something*)(n->value);
    Mesh_viewer_faces *mf=(Mesh_viewer_faces*)(ms->evolution);
    mf->color_rows=mesh.num_v(&mesh);
    mf->random_color(mf);
    //mf->normal_rows=mesh.num
    get_data_from_3dim_cell(&mesh,&(mf->Data),&(mf->Data_index),&(mf->Data_rows),&(mf->Data_index_rows));
    mf->normal_rows=mf->Data_rows;
*/
    /*char sp_name[]="faces";
    Node* n=mw.create_something(&mw,sp_name);
    auto ms=(Mesh_viewer_something*)(n->value);
    auto mf=(Mesh_viewer_faces*)(ms->evolution);
    mf->color_rows=mesh.num_c(&mesh);
    mf->random_color(mf);
//    get_data_from_2dim_cell(&mesh,&(mf->Data),&(mf->Data_index));
    mf->normal_rows=mesh.num_c(&mesh);
    mf->Data_rows=mesh.num_v(&mesh);
    mf->Data_index_rows=mesh.num_c(&mesh);
    free_node(n);*/
    //test_camera_and_intera(&mw);
    //Mesh_viewer_opengl_Interpreter moi;
    //Mesh_viewer_opengl_Interpreter_init(&moi);
    //moi.world=&mw;
    //moi.routine_show(&moi);
    
}

void once(Mesh_viewer_Intera*mi)
{
    myani*anima=(myani*)(mi->representation);
    Mesh* mesh=anima->mesh;

    Mesh_viewer_faces* mf=anima->mf;
    anima->i=mesh->cell_id;
    if(mi->g_info->key_action==MESH_VIEWER_PRESS&&mi->g_info->key==MESH_VIEWER_KEY_A)
    {
        Node* node_f=NULL;
        for(auto f_it=mesh->faces.begin();f_it!=mesh->faces.end();f_it++)
        {  
            int *f_id=(int*)malloc(sizeof(int));
            if(!libcell_face_is_convex(mesh,f_it->second))
            {
                *f_id=f_it->second->id;
                node_f=node_overlying(node_f,f_id);
            }
        }

        printf("node_f:%d\n",node_size(node_f));
        template_f* f0=NULL;
        Node* n_it=node_f;       
        auto iter=mesh->faces.find(*((int*)(n_it->value)));
        if(iter!=mesh->faces.end())
        {
            f0=iter->second;
            bool flag=convex_adjust(mesh,f0);
        }
        free_node_value(node_f);
        free_node(node_f); 
        Node* node_c=NULL;
        for(auto c_it=mesh->cells.begin();c_it!=mesh->cells.end();c_it++)
        {
            if(c_it->second->id>=anima->i)
            {
                node_c=node_overlying(node_c,c_it->second);
            }
        }
        n_it=node_c;
        int f_size=node_size(node_c);
        mf->Data_rows=f_size*3;
        mf->Data_index_rows=f_size;
        mf->color_rows=mf->Data_index_rows;
        float color[3]={1.0,0.0,0.0};
        mf->set_color(mf,color);
        mf->normal_rows=mf->Data_index_rows;
        mf->Data=(float*)malloc(sizeof(float)*3*mf->Data_rows);
        mf->Data_index=(unsigned int*)malloc(sizeof(unsigned int)*mf->Data_index_rows*4);
        template_c* c0=NULL;int i=0,k=0;
        while(n_it!=NULL)
        {
            c0=(template_c*)(n_it->value);
            mf->Data_index[i]=3;
            i++;
            for(auto iter=mesh->cv_begin(mesh,*c0);iter!=mesh->cv_end(mesh,*c0);iter++)
            {
                mf->Data_index[i]=k;
                for(int j=0;j<3;j++)
                {
                    mf->Data[j+k*3]=(*iter).point[j];
                }
                i++;k++; 
            }
            n_it=(Node*)(n_it->Next);
        }

        free_node(node_c);

        anima->moi->update_data(anima->moi); 


    }
/*   if(mi->g_info->key_action==MESH_VIEWER_PRESS&&mi->g_info->key==MESH_VIEWER_KEY_A)
   {
        printf("a\n");
        template_v* v1=mesh->create_vertexv(mesh,(anima->v)[anima->i],3);
        increasing_convex_hull(mesh,v1);
        anima->i++;
        mesh->printself(mesh);
        //update
        
        free(mf->Data);
        mf->Data=0;
        free(mf->Data_index);
        mf->Data_index=0;
        mf->color_rows=mesh->num_v(mesh);
        mf->random_color(mf);
        get_data_from_2dim_cell(mesh,&(mf->Data),&(mf->Data_index));
        mf->normal_rows=mesh->num_c(mesh);
        mf->compute_normal(mf);
        mf->Data_rows=mesh->num_v(mesh);
        mf->Data_index_rows=mesh->num_c(mesh);
        anima->moi->update_data(anima->moi);  
   }
   else if(mi->g_info->key_action==MESH_VIEWER_PRESS&&mi->g_info->key==MESH_VIEWER_KEY_B)
   {
   }*/



}
void test_camera_and_intera(Mesh_viewer_world* mw)
{
    char camera[]="Camera";

    Node*n=mw->create_something(mw,camera);

    auto ms=(Mesh_viewer_something*)(n->value);

    Mesh_viewer_camera* mc=(Mesh_viewer_camera*)(ms->evolution); mc->is_using=1;

    Matrix4x4* p=Projection(M_PI/3.0f,(float)(mw->g_info->resolution[0])/(float)(mw->g_info->resolution[1]),0.5f,200.0f);
   // p->print_self(p);
    Matrix4x4_copy_data_float(mc->Proj,p);
    Matrix4x4_free(p);    
    free_node(n);
    n=0;

    char Intera[]="Intera";
    n=mw->create_something(mw,Intera);
    ms=(Mesh_viewer_something*)(n->value);
    Mesh_viewer_Intera* mi=(Mesh_viewer_Intera*)(ms->evolution);
    //mi->g_info=mw->g_info;
    //mi->prop=(void*)mc;
    Mesh_viewer_Arcroll* ma=(Mesh_viewer_Arcroll*)malloc(sizeof(Mesh_viewer_Arcroll));
    Mesh_viewer_Arcroll_init(ma);
    mi->representation=(void*)ma;
    ma->mc=mc;
    mi->cursor_position_callback=Mesh_viewer_Arcroll_cursor_position_callback;
    mi->scroll_callback=Mesh_viewer_Arcroll_scroll_callback;
    free_node(n);

}

void test_convex()
{

    srand((unsigned)time(0));
    Mesh mesh;
    A_Mesh_init(&mesh);
    double **v=(double**)malloc(sizeof(double*)*100);
    for(int i=0;i<100;i++)
    {
        v[i]=(double*)malloc(sizeof(double)*3);

    }
/*    v[0][0]=-0.5;v[0][1]=-0.5;v[0][2]=0.5;
    v[1][0]=0.5;v[1][1]=-0.5;v[1][2]=0.5;
    v[2][0]=-0.5;v[2][1]=0.5;v[2][2]=0.5;
    v[3][0]=-0.5;v[3][1]=0.5;v[3][2]=-0.5;
    v[4][0]=0.5;v[4][1]=0.5;v[4][2]=0.5;
    v[5][0]=0.5;v[5][1]=0.5;v[5][2]=-0.5;
    v[6][0]=-0.5;v[6][1]=-0.5;v[6][2]=-0.5;
    v[7][0]=0.5;v[7][1]=-0.5;v[7][2]=-0.5;*/
    v[0][0]=-0.026015;v[0][1]= 0.112578;v[0][2]= 0.036387; 
    v[1][0]=-0.032178;v[1][1]=0.174119;v[1][2]= -0.002633; 
    v[2][0]=-0.080718;v[2][1]=0.152855;v[2][2]= 0.030245; 
    v[3][0]=-0.023129;v[3][1]= 0.112186;v[3][2]= 0.038644; 
    v[4][0]=0.016493;v[4][1]= 0.122293;v[4][2]= 0.031423; 
    v[5][0]=-0.024837;v[5][1]= 0.156574;v[5][2]= -0.003775; 
    v[6][0]=0.045263;v[6][1]= 0.086356;v[6][2]= 0.020037; 
    v[7][0]=-0.071339;v[7][1]= 0.155715;v[7][2]= 0.007126; 
    v[8][0]=-0.010076;v[8][1]= 0.125949;v[8][2]= 0.029738; 
    v[9][0]=-0.072013;v[9][1]= 0.154518;v[9][2]= 0.030398; 
    v[10][0]=-0.074228;v[10][1]= 0.163163;v[10][2]= -0.014495; 
    v[11][0]=-0.075178; v[11][1]=0.166521;v[11][2]= -0.022380; 
    v[12][0]=-0.077120;v[12][1]= 0.165160;v[12][2]= -0.023763; 
    v[13][0]=-0.015844;v[13][1]= 0.114601;v[13][2]= 0.038775; 
    v[14][0]=-0.076832;v[14][1]= 0.134852;v[14][2]= 0.051454; 
    v[15][0]=-0.017349;v[15][1]= 0.110376;v[15][2]= 0.041041; 
    v[16][0]=-0.036221;v[16][1]= 0.170661;v[16][2]= 0.000806; 
    v[17][0]=-0.037684;v[17][1]= 0.171196;v[17][2]= 0.000807; 
    v[18][0]=-0.091535;v[18][2]= 0.145199;v[18][2]= 0.023217; 
    v[19][0]=0.008183;v[19][1]= 0.131521;v[19][2]= 0.014497; 
    v[20][0]=-0.057447;v[20][1]= 0.041722;v[20][2]= 0.046094; 
    v[21][0]=-0.052683;v[21][1]= 0.043904;v[21][2]= 0.045084; 
    v[22][0]=-0.064144;v[22][1]= 0.161751;v[22][2]= -0.019570; 
    v[23][0]=-0.030386;v[23][1]= 0.178785;v[23][2]= -0.005728; 
    v[24][0]=-0.025781;v[24][1]= 0.115463;v[24][0]= 0.033662; 
    v[25][0]=0.003016;v[25][1]= 0.042381;v[25][2]= 0.045895; 
    v[26][0]=-0.024994;v[26][1]= 0.180249;v[26][2]= -0.009198; 
    v[27][0]=-0.059546;v[27][1]= 0.156361;v[27][2]= 0.007485;
    v[28][0]=-0.062224;v[28][1]=0.155123;v[28][2]= 0.005565; 
    convex_init(&mesh,v,3);
    //   from_v_createconvex(&mesh,v,20,3);
    for(auto iter=mesh.cells.begin();iter!=mesh.cells.end();iter++)
    {
        for(auto it_cv=mesh.cv_begin(&mesh,*(iter->second));it_cv!=mesh.cv_end(&mesh,*(iter->second));it_cv++)
        {
            //printf("id: %d  ",(*it_cv).id);
        }
    }
    Mesh_viewer_world mw;
    Mesh_viewer_world_init(&mw);
    char sp_name[]="faces";
    Node* n=mw.create_something(&mw,sp_name);
    auto ms=(Mesh_viewer_something*)(n->value);
    auto mf=(Mesh_viewer_faces*)(ms->evolution);
    mf->color_rows=mesh.num_v(&mesh);
    mf->random_color(mf);
    //float temp_color[]={0.5,0.5,0.1};
    //mf->set_color(mf,temp_color);
    get_data_from_2dim_cell(&mesh,&(mf->Data),&(mf->Data_index));
    mf->normal_rows=mesh.num_c(&mesh);
    mf->Data_rows=mesh.num_v(&mesh);
    mf->Data_index_rows=mesh.num_c(&mesh);
    free_node(n);
    test_camera_and_intera(&mw);
    char intera[]="Intera";
    n=mw.create_something(&mw,intera);
    ms=(Mesh_viewer_something*)(n->value);
    auto mi=(Mesh_viewer_Intera*)(ms->evolution);
    myani anima;
    anima.mesh=&mesh;anima.mf=mf;anima.i=4;anima.v=v;
    //mi->g_info=mw.g_info;
    mi->representation=(void*)(&anima);

    mi->key_callback=once;
    free_node(n);
    Mesh_viewer_opengl_Interpreter moi;
    Mesh_viewer_opengl_Interpreter_init(&moi);
    anima.moi=&moi;
    moi.world=&mw;
    moi.routine_show(&moi);
    

}
static void examine(Mesh* mesh)
{
    Node* node=NULL;
    for(auto iter=mesh->faces.begin();iter!=mesh->faces.end();iter++)
    {
        if(mesh->face_is_boundary(mesh,iter->second))
        {
            template_c* c=0;
            if(iter->second->halffaces[0]->cell!=0)
            {
                c=(template_c*)(iter->second->halffaces[0]->cell);
            }
            else if(iter->second->halffaces[1]->cell!=0)
            {
                c=(template_c*)(iter->second->halffaces[1]->cell);
            }
            if(c!=0)
            {

                node=node_overlying(node,(void*)c);
            }
        
        } 
    }
    while(node!=0)
    {
       mesh->delete_cell(mesh,*((template_c*)(node->value)),true);
       node=(Node*)(node->Next); 
    }

}
void test3D(){
    Mesh mesh,mesh1;
    A_Mesh_init(&mesh);A_Mesh_init(&mesh1);
    _ReadCell_(&mesh,"delauny_sphere.cell");
    //examine(&mesh);
    //_ReadCell_(&mesh1,"face.cell");
    _ReadCell_(&mesh1,"delauny_circle.cell");
/*    double temp_sum=0;
    for(auto iter=mesh1.vertices.begin();iter!=mesh1.vertices.end();iter++)
    {
        temp_sum=0;
        double* temp_v=(double*)malloc(sizeof(double)*(iter->second->point_size+1));
        for(int j=0;j<iter->second->point_size;j++)
        {
            temp_v[j]=iter->second->point[j];
            temp_sum+=iter->second->point[j]*iter->second->point[j];
        }
        temp_v[iter->second->point_size]=temp_sum/2.0;
        iter->second->point_size++;
        free(iter->second->point);
        iter->second->point=temp_v;
    }
*/
    //mesh1.printself(&mesh1);
        // _ReadOff_(&mesh,"cube.off",3);
    printf("%d %d %d\n",mesh.num_v(&mesh),mesh.num_f(&mesh),mesh.num_c(&mesh));
    Mesh_viewer_world mw;
    Mesh_viewer_world_init(&mw);
    char ch[]="faces,faces";
    Node* n=mw.create_something(&mw,ch);
    Mesh_viewer_something*ms=(Mesh_viewer_something*)(n->value);
    Mesh_viewer_faces *mf=(Mesh_viewer_faces*)(ms->evolution);
    mf->color_rows=mesh.num_v(&mesh);
    //mf->is_reversal_normal=1; 
    mf->random_color(mf);
    get_data_from_3dim_cell(&mesh,&(mf->Data),&(mf->Data_index),&(mf->Data_rows),&(mf->Data_index_rows)); 
    mf->normal_rows=mf->Data_rows;
    ms->disappear=1;
    n=(Node*)(n->Prev);
    ms=(Mesh_viewer_something*)(n->value);
    mf=(Mesh_viewer_faces*)(ms->evolution);
    mf->color_rows=mesh1.num_c(&mesh1);
    mf->random_color(mf);
    //float temp_color[]={0.1,0.9,0.3};
    //mf->set_color(mf,temp_color);
    mf->normal_rows=mesh1.num_c(&mesh1);
    get_data_from_2dim_cell(&mesh1,&(mf->Data),&(mf->Data_index));
    mf->Data_rows=mesh1.num_v(&mesh1);
    mf->Data_index_rows=mesh1.num_c(&mesh1);
   // ms->disappear=1;
        //get_array_from_libcell_mesh(&mesh,&(mf->Data),&(mf->Data_index)); 
        //mf->Data_rows=mesh.num_v(&mesh);
        //mf->Data_index_rows=mesh.num_c(&mesh);
        //mf->normal_rows=mesh.num_c(&mesh);

    free_node(n);
    char edges[]="edges";
    n=mw.create_something(&mw,edges);
    ms=(Mesh_viewer_something*)(n->value);
    Mesh_viewer_edges*me=(Mesh_viewer_edges*)(ms->evolution);
    me->Data_rows=4;
    me->Data_index_rows=2;

    float*v=(float*)malloc(sizeof(float)*3*me->Data_rows);
    unsigned int*f=(unsigned int*)malloc(sizeof(unsigned int)*2*me->Data_index_rows);
    v[0]=1.0;v[1]=0.0;v[2]=0.0;
    v[3]=1.0;v[4]=-0.0;v[5]=0.5;
     v[6]=-1.0;v[7]=0.0;v[8]=0.0;
    v[9]=-1.0;v[10]=-0.0;v[11]=0.5;

    f[0]=0;f[1]=1;f[2]=2;f[3]=3;
    me->Data=v;
    me->Data_index=f;
    me->color_rows=4;
    me->color=(float*)malloc(sizeof(float)*3*me->color_rows);
    //me->color[0]=1.0;me->color[1]=0.0;me->color[2]=0.0;
    //me->color[3]=1.0;me->color[4]=0.0;me->color[5]=0.0;

    free_node(n);
    test_camera_and_intera(&mw);
        //mf->compute_normal(mf);
    Mesh_viewer_opengl_Interpreter moi;
    Mesh_viewer_opengl_Interpreter_init(&moi);
    moi.world=&mw;
    moi.routine_show(&moi);
    Mesh_free(&mesh);
//    moi.render(&moi);

}
void convex_init(Mesh* mesh,double **VV,int cols)
{

     mesh->simplex=1;
    for(int i=0;i<cols+1;i++)
    {
        mesh->create_vertexv(mesh,VV[i],cols);

    }
    int *index=(int*)malloc(sizeof(int)*(cols+1));
    for(int i=0;i<cols+1;i++)
    {
        index[i]=i;
    }
    if(ori_area_simplex(VV,cols+1,cols)<0)
    {
        index[0]=1;
        index[1]=0;
   
    }
  template_v **temp_v=(template_v**)malloc(sizeof(template_v*)*cols);

    for(int i=0;i<cols+1;i++)
    {
        int k=0;
//        template_v* temp_v[cols];
        for(int j=0;j<cols+1;j++)
        {
            if(j!=i)
            {
                temp_v[k]=mesh->vertices.find(index[j])->second;
                k++;
            }

        }
        if((cols-i)%2!=0)
        {
            template_v* temp_v1=temp_v[0];
            temp_v[0]=temp_v[1];
            temp_v[1]=temp_v1;
        }
        mesh->create_cellv(mesh,temp_v,cols);
        for(int j=0;j<cols;j++)
        {
            printf("%d ",(*(temp_v[j])).id); 
        }
        printf("\n");
    }
    free(temp_v);

    free(index);
}
