#include "include/opt.h"

#include "include/Cholesky.h"
#include<iostream>
void Libcell_Opt_init(L_Opt* op)
{
    op->node_v=NULL;
    op->v_size=0;
    op->v0=NULL;
    op->mesh=NULL;
    op->sum_meas=0;
    op->H=NULL;
    op->Wi=NULL;
    //Ai不变
    op->Ai=NULL;
    op->rows=0;
}
//用到了c的prop属性,没有用到v的prop
void libcell_compute_dual_point(template_m*mesh,template_c* c_)
{
    int dim=mesh->dimension;
    Eigen::MatrixXd A,temp_re,temp_b;
    temp_b.resize(dim,1);
    A.resize(dim,dim);
    template_v*v0,*vi;
  
    double *temp_v=(double*)malloc(sizeof(double)*(dim));
    c_->prop=(void*)temp_v;
    v0=(template_v*)(c_->vertices[0]);
    for(int i=1;i<c_->vertices_size;i++)
    {
        vi=(template_v*)(c_->vertices[i]);
        for(int j=0;j<vi->point_size-1;j++)
        {
                A.coeffRef(i-1,j)=vi->point[j]-v0->point[j];
            
        }
           temp_b.coeffRef(i-1,0)=vi->point[dim]-v0->point[dim]; 
    }

    temp_re=A.inverse()*temp_b; 
    for(int i=0;i<dim;i++)
    {
        temp_v[i]=temp_re.coeff(i,0);
    }
   
}
static Node* two_points_cells(template_m*m,template_v*v0,template_v* v1)
{
    Node* re=NULL;
    Node* l=v1->cells;
    for(auto iter=m->vc_begin(m,*v0);iter!=m->vc_end(m,*v0);iter++)
    {
        if(node_find(l,quote(iter))!=NULL)
        {
            re=node_overlying(re,quote(iter)); 
        } 
    }
    return re;
}
static void test_cut_num_(double **v,int rows,int cols)
{
    for(int i=0;i<rows;i++)
    {
        for(int j=0;j<cols;j++)
        {
            int a=1000000*v[i][j];
            v[i][j]=a/1000000.0; 
        }
    
    }
}

static double two_points_dual_data(template_m*m,template_v*v0,template_v*v1)
{
    double re=0;
    Node* n_it=two_points_cells(m,v0,v1);
    Node* n_it1=n_it;
    int dim=m->dimension;
    int cols=dim;
    int rows=node_size(n_it);
    
    double **v=(double**)malloc(sizeof(double*)*rows);
    for(int i=0;i<rows;i++)
    {
        v[i]=(double*)malloc(sizeof(double)*cols);
    }
    template_c*c0=NULL;double *p=NULL;
    int k=0;
    while(n_it!=NULL)
    {
        c0=(template_c*)(n_it->value);
        p=(double*)(c0->prop);
        for(int j=0;j<cols;j++)
        {
            v[k][j]=p[j]; 
        }
        k++;
        n_it=(Node*)(n_it->Next);
    }
    /*double re1=0;
    for(int j=0;j<cols;j++)
    {
         re1+=(v[0][j]-v[1][j])*(v[0][j]-v[1][j]); 
    }*/
    test_cut_num_(v,rows,cols);
    re=compute_convex_area(v,rows,cols,dim-1);
    //printf("bijiao:%lf   %lf\n",re,sqrt(re1));
    double norm=0;
    for(int j=0;j<cols;j++)
    {
        norm+=(v0->point[j]-v1->point[j])*(v0->point[j]-v1->point[j]);
    }
    re=re/(sqrt(norm));
    for(int i=0;i<rows;i++)
    {
        free(v[i]);
    }
    free(v);
    free_node(n_it1);
    return re;
}

static double get_measure(L_Opt* myopt)
{    Mesh* mesh=myopt->mesh;
    template_v*v0=myopt->v0,*v1=NULL;
    double re=0;
    int dim=mesh->dimension;
    int cols=dim+1,rows=0;
    Node* n_it=NULL,*n_it1=NULL;
    double*p0=NULL;
    template_c* c0=NULL;
    double temp_v1;
    n_it1=myopt->node_v;
    while(n_it1!=NULL)
    {
        v1=(template_v*)(n_it1->value);

        n_it=v1->cells;
        rows=node_size(n_it);
        double **v=(double**)malloc(sizeof(double*)*rows);
        for(int i=0;i<rows;i++)
        {
            v[i]=(double*)malloc(sizeof(double)*cols);
            
        }
        int i=0;
        while(n_it!=NULL)
        {
            c0=(template_c*)(n_it->value);
            p0=(double*)(c0->prop);
            for(int j=0;j<dim;j++)
            {
                v[i][j]=p0[j];
            }
            v[i][dim]=0;
            n_it=(Node*)(n_it->Next);
            i++;
        }
        //here
	/*	if(rows>=15)
		{
            
            test_cut_num_(v,rows,cols);
		    printf("h rows:%d\n",rows);
            for(int i=0;i<rows;i++)
            {
                for(int j=0;j<cols;j++)
                {
                    printf("%.30f,",v[i][j]);
                }
                printf("\n");
            }
        }
		else
		{
		printf("rows:%d\n",rows);
            for(int i=0;i<rows;i++)
            {
                for(int j=0;j<cols;j++)
                {
                    printf("%.30f,",v[i][j]);
                }
                printf("\n");
            }
		}*/	    
	    test_cut_num_(v,rows,cols);
        temp_v1=compute_convex_area(v,rows,cols,dim);
        myopt->Wi[*((int*)(v1->prop1))]=temp_v1;
        printf("wi:%lf \n",temp_v1);
        re+=temp_v1;
        for(i=0;i<rows;i++)
        {
            free(v[i]);
        }
        free(v);

        n_it1=(Node*)(n_it1->Next);
    }
    myopt->sum_meas=re;
    return re;
}
static double compute_energy(L_Opt*myopt)
{
//    int rows=myopt->Wi.rows();
    double re=0;
    for(int i=0;i<myopt->rows;i++)
    {
        printf("wi:%lf\n",myopt->Wi[i]);
        re+=myopt->Wi[i]*myopt->Wi[i];
    }
    printf("energy:%lf\n",re);
    return re;

}
static void fill_node_v(L_Opt*myopt)
{
    Mesh* mesh=myopt->mesh;
    template_v*v0=myopt->v0;
    if(myopt->node_v!=NULL)
    {
        free_node(myopt->node_v);
        myopt->node_v=NULL;
    }

    for(auto iter=mesh->vertices.begin();iter!=mesh->vertices.end();iter++)
    {
        if(iter->second==v0)
        {
            continue;
        }
        if(iter->second->cells==0)
        {
            continue; 
        }
        myopt->node_v=node_overlying(myopt->node_v,iter->second);
    } 

}
static void sort_node_v(L_Opt*myopt)
{
    template_v* v=NULL;
    Node* n_it=myopt->node_v;
    int i=0;
    int*id=NULL;
    while(n_it!=NULL)
    {
        v=(template_v*)(n_it->value);
        id=(int*)(v->prop1);
        *id=i;
        n_it=(Node*)(n_it->Next);
        i++;
    }
    printf("sort:%d\n",i);

}
static void fill_H(L_Opt* myopt)
{
     
    template_v*v=NULL,*v1=NULL,*v0=NULL;
    v0=myopt->v0;
    Node* n_it=myopt->node_v;
    Mesh*m=myopt->mesh;
    double**H=myopt->H;
    int i=0,j=0;
    double temp_sum=0;
  
    while(n_it!=NULL)
    {
    
        v=(template_v*)(n_it->value); 
        Node* iter=m->vv_begin(m,*v),*n_it1=NULL;
        n_it1=iter;
        i=*((int*)(v->prop1));
        temp_sum=0;
        while(iter!=NULL)
        {
            v1=(template_v*)(iter->value);
            if(v1!=v0)
            {

                j=*((int*)(v1->prop1));
                if(H[j][i]!=0)
                {
                    H[i][j]=H[j][i];
                }
                else 
                {
                    H[i][j]=two_points_dual_data(m,v,v1);
                }
                temp_sum-=H[i][j];
            }
            iter=(Node*)(iter->Next);
        }
        H[i][i]=temp_sum;

        free_node(n_it1);
        n_it=(Node*)(n_it->Next);
    }

}
void iteration(L_Opt* myopt)
{
    fill_node_v(myopt);
    sort_node_v(myopt);
    int rows=node_size(myopt->node_v);
    myopt->rows=rows;
    double* Wi=(double*)malloc(sizeof(double)*rows);
    myopt->Wi=Wi;
        //fill Wi

    get_measure(myopt);
    printf("sum_meas:%lf\n",myopt->sum_meas);
    Node* n_it=myopt->node_v;
    while(n_it!=NULL)
    {
        template_v*v=(template_v*)(n_it->value);
        int i=*((int*)(v->prop1));
        Wi[i]/=myopt->sum_meas;
        Wi[i]=myopt->Ai[v->id]-Wi[i];
        n_it=(Node*)(n_it->Next); 
    }

     printf("rows:%d meas:%lf\n",rows,myopt->sum_meas);
    double**H=(double**)malloc(sizeof(double*)*rows);
    for(int i=0;i<rows;i++)
    {
        H[i]=(double*)malloc(sizeof(double)*rows);
        memset(H[i],0,sizeof(double)*rows);
    }
    myopt->H=H; 
    double **L=(double**)malloc(sizeof(double*)*rows);
    for(int i=0;i<rows;i++)
    {
        L[i]=(double*)malloc(sizeof(double)*rows);
        memset(L[i],0,sizeof(double)*rows);
    }
    for(int i=0;i<rows;i++)
    {
        L[i][i]=1;
    }
    double* D=(double*)malloc(sizeof(double)*rows);

    fill_H(myopt);
    cholesky_decomp(H,L,D,rows-1);
    double*x=cholesky_solve(H,L,D,Wi,rows-1);    
    n_it=myopt->node_v;
    int i=0;double t=0.014;
    while(n_it!=NULL)
    {
        template_v*v=(template_v*)(n_it->value);
        v->point[myopt->mesh->dimension]+=x[i]*t;
        n_it=(Node*)(n_it->Next); 
        i++;
    }
    compute_energy(myopt);
    
}

/*static double data_from_two_halffaces(template_m*m,template_hf*hf1,template_hf*hf2,bool&is_linked)
{
   template_v* v0=NULL;
   int k=0;
   int flag=1;
   for(int i=0;i<hf1->vertices_size;i++)
   {
       flag=1;
        for(int j=0;j<hf2->vertices_size;j++)
        {
            if(hf1->vertices[i]==hf2->vertices[j])
            {
                flag=0;
                break;
            }
        
        }
        if(flag==1)
        {
            k++;
            v0=(template_v*)(hf1->vertices[i]);
        }
   }
    if(k==1)
    {
        is_linked=true;
    
    }
    else
    {
        is_linked=false;
        return 0;
    }
    int cols=v0->point_size;
    int rows=cols+1;
    double **VV=(double**)malloc(sizeof(double*)*rows);
    for(int i=0;i<rows;i++)
    {
        VV[i]=(double*)malloc(sizeof(double)*cols);
    }
    for(int j=0;j<cols;j++)
    {

        VV[cols][j]=v0->point[j];
    }
    for(int i=0;i<cols;i++)
    {
        for(int j=0;j<cols;j++)
        {
            VV[i][j]=((template_v*)(hf2->vertices[i]))->point[j]; 
        }
    }
    return ori_area_simplex(VV,rows,cols);

}
static Node* mesh_is_not_convex_halfface(template_m*m,Node* node_hf)
{
//   m->external_cell_init_(m);
    Node* re=NULL;
    template_hf *hf1=NULL,*hf2=NULL;
    bool is_linked;
    double re1;
   while(node_hf!=NULL)
   {
       Node* it=(Node*)(node_hf->Next);
       while(it!=NULL)
       {
           hf1=(template_hf*)(node_hf->value);
           hf2=(template_hf*)(it->value);
            //Node* temp=m->intersection_two_faces(m,hf1->face,hf2->face);
           re1=data_from_two_halffaces(m,hf1,hf2,is_linked); 
           if(is_linked&&re1<0)
            {
              re=node_overlying(re,hf1);
              //re=node_overlying(re,hf2);
              return re;
            }
            it=(Node*)(it->Next);
       }
    node_hf=(Node*)(node_hf->Next);
   }
    return re;
}*/
static bool halfface_is_special(template_m*m,template_hf*hf)
{
    template_hf* hf1=m->s_opposite_halfface(hf);
    if(hf1->cell==NULL)
    {
        return true;
    }
    else    
    {
        if(((template_c*)(hf1->cell))->prop1==NULL&&((template_c*)(hf->cell))->prop1!=NULL)
        {
            return true;      
        }
        /*else if()
        {
        }*/
        
    }
    return false;
}
static Node* extract_halfface_from_cell(template_m*m,Node*node_c)
{
    printf("begin extract half\n");

    Node* node_hf=NULL;Node* it_c=node_c;
    while(it_c!=NULL)
    {
    
        for(auto iter=m->chf_begin(m,*((template_c*)(it_c->value)));iter!=m->chf_end(m,*((template_c*)(it_c->value)));iter++)
        {
            if(halfface_is_special(m,quote(iter)))
            {
                for(int j=0;j<quote(iter)->vertices_size;j++)
                {
                    printf(" id :%d ",((template_v*)(quote(iter)->vertices[j]))->id);
                }
                printf("\n");
                node_hf=node_overlying(node_hf,quote(iter));
            }
        
        }
        it_c=(Node*)(it_c->Next);
    }
    printf("\n end extract half\n");

    return node_hf;
}
static Node* extract_vertices_from_cell(template_m*m,Node* node_c)
{
    printf("begin extract vert\n");
    Node* re=NULL;
    Node* it_c=node_c;
    template_c*c=NULL;
    template_v*v0=NULL;
    while(it_c!=NULL)
    {
        c=(template_c*)(it_c->value);
        for(auto iter=m->cv_begin(m,*c);iter!=m->cv_end(m,*c);iter++)
        {
            v0=quote(iter);

            if(node_find(re,v0)==NULL)
            {
                printf("v id :%d  ",v0->id);
                re=node_overlying(re,v0);
            }
        }
        it_c=(Node*)(it_c->Next);
    }
    
    printf("\n end extract vert\n");

    return re;
}
static bool jianyan_faces(template_m*m)
{
    printf("jiance\n");
    for(auto iter=m->faces.begin();iter!=m->faces.end();iter++)
    {
        if(iter->second==NULL)
        {
            printf("chuwenti f==NULL\n");
            return false;
        }
        if(m->face_is_boundary(m,(iter->second)))
        {
            
            printf("chuxian boundary id:%d\n",iter->second->id);
            for(int j=0;j<iter->second->vertices_size;j++)
            {
                printf("%d ",((template_v*)(iter->second->vertices[j]))->id); 
            }
            printf("\n");
            return false;
        }
    }

    return true;
}
static Node*  get_not_convex_hf(template_m*m,Node* node_c)
{
    printf("begin get not convex hf\n");
    Node* node_hf=extract_halfface_from_cell(m,node_c);
    Node* node_v=extract_vertices_from_cell(m, node_c);

    Mesh mesh;
    A_Mesh_init(&mesh);
    Node*n_it=node_v; template_v*v0=NULL,*v1=NULL,*v2=NULL;
    int cols=m->dimension+1;
    double* temp_v=(double*)malloc(sizeof(double)*cols);
    memset(temp_v,0,sizeof(double)*cols);
    temp_v[cols-1]=10000.0;
    v2=mesh.create_vertexv(&mesh,temp_v,cols);
    while(n_it!=NULL)
    {
        v1=(template_v*)(n_it->value);
        for(int j=0;j<cols;j++)
        {
            temp_v[j]=v1->point[j]; 
        }
        v0=mesh.create_vertexv(&mesh,temp_v,cols);
        v0->prop=(void*)v1;
        n_it=(Node*)(n_it->Next);
    }
    mesh_createconvex(&mesh);
    mesh.delete_vertex(&mesh,*v2,true);
    mesh.external_cell_init_(&mesh);
    Node* c_it=mesh.external_cell.halffaces;
    //template_c* c0=NULL;
    template_hf* hf0=NULL;template_f* f0=NULL;
    template_v**temp_v1=(template_v**)malloc(sizeof(template_v*)*(m->dimension));

    printf("test half\n");
    while(c_it!=NULL)
    {   
       hf0=(template_hf*)(c_it->value);
       for(int j=0;j<m->dimension;j++)
       {
        temp_v1[j]=(template_v*)(((template_v*)(hf0->vertices[j]))->prop);
        printf("vid :%d  ",temp_v1[j]->id);
       }
       printf("\n");

        f0=m->get_facev(m,temp_v1,m->dimension);
        if(f0!=NULL)
        {
            f0->halffaces[0];
            node_hf=node_delete_value(node_hf,f0->halffaces[0]);
            node_hf=node_delete_value(node_hf,f0->halffaces[1]); 
        }
        c_it=(Node*)(c_it->Next);
    }

    printf("test half\n");
    template_v**temp_v2=(template_v**)malloc(sizeof(template_v*)*(m->dimension+1));
    c_it=node_c;

    jianyan_faces(m);
    if(node_hf==NULL)
    {
        printf("begin tianchong\n");
        printf("delete cell v id\n");
        while(c_it!=NULL)
        {
            for(int j=0;j<((template_c*)(c_it->value))->vertices_size;j++)
            {
                printf("%d ",((template_v*)(((template_c*)(c_it->value))->vertices[j]))->id);
            
            }
            printf("\n");
            m->delete_cell(m,*((template_c*)(c_it->value)),true); 
            c_it=(Node*)(c_it->Next);
        }

        printf("create cell v id\n");
        for(auto iter1=mesh.cells.begin();iter1!=mesh.cells.end();iter1++)
        {
            for(int j=0;j<iter1->second->vertices_size;j++)
            {
                temp_v2[j]=(template_v*)(((template_v*)(iter1->second->vertices[j]))->prop);
                printf("%d ",temp_v2[j]->id);
            }

            template_c*c0=m->create_cellv(m,temp_v2,m->dimension+1);
            printf("\n cell id:%d\n",c0->id);
        }
       // printf("create dizhi:%p\n",m->prop);
        printf("end tianchong\n");

    }

    jianyan_faces(m);
    Mesh_free(&mesh);
    free(temp_v);
    free(temp_v1);
    free(temp_v2);
    free_node(node_v);
    printf("end get not convex hf\n");

    return node_hf;
}
//检验vv_begin
bool convex_adjust(template_m*m,template_f* f)
{
    printf("begin convex adjust\n");
   // jianyan_faces(m);

    if(f==NULL||m==NULL)
    {
        return false;
    }
    //f不能是边界
    if(m->simplex!=1)
    {
        printf("m is not simplex\n");
        return false;
    }
    if(m->face_is_boundary(m,f))
    {
        return false;  
    }
    int * flag=(int*)malloc(sizeof(int));
    *flag=1;
    Node* node_c=NULL;
    ((template_c*)(f->halffaces[0]->cell))->prop1=(void*)flag;
    node_c=node_overlying(node_c,f->halffaces[0]->cell);
    ((template_c*)(f->halffaces[1]->cell))->prop1=(void*)flag;
    node_c=node_overlying(node_c,f->halffaces[1]->cell);
     
    Node* node_hf=get_not_convex_hf(m,node_c);
    printf("my\n");
    jianyan_faces(m);
    printf("my\n");

    template_c* c0=NULL;template_hf* hf0=NULL;
    while(node_hf!=NULL)
    {
        printf("xunhuan\n");
        jianyan_faces(m);
        printf("xunhuan here1\n");
        hf0= ((template_hf*)(node_hf->value));

        hf0=m->s_opposite_halfface(hf0);

        c0=(template_c*)(hf0->cell);

        printf("xunhuan here2\n");
        if(c0==NULL)
        {
            printf("here error:\n");
        }
 
        if(c0->prop1!=NULL)
        {
            printf("here error\n");
        }
        printf("xunhuan here3\n");

        c0->prop1=(void*)flag;
        node_c=node_overlying(node_c,c0);
        free_node(node_hf);
        printf("xunhuan here4\n");
        node_hf=get_not_convex_hf(m,node_c); 
    }
    //printf("here\n");
    //jianyan_faces(m);
    //jianyan_faces(m);

    free_node(node_c);
    printf("here1\n");
     jianyan_faces(m);

    free(flag);
    printf("end convex adjust\n");

   return true;
}


