#ifndef MESH_VIEWER_ARCROLL
#define MESH_VIEWER_ARCROLL
#include<Mesh_Viewer_World_Manager/Mesh_viewer_manager.h>
//#include<Mesh_view/Mesh_viewer_interpreter.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct Mesh_viewer_Arcroll
{
    float* old_mouse_coord;
	Mesh_viewer_camera*mc;
}Mesh_viewer_Arcroll;
void Mesh_viewer_Arcroll_init(Mesh_viewer_Arcroll*);
void Mesh_viewer_Arcroll_cursor_position_callback(Mesh_viewer_Intera*);
void Mesh_viewer_Arcroll_mouse_button_callback(Mesh_viewer_Intera*);
void Mesh_viewer_Arcroll_framebuffer_size_callback(Mesh_viewer_Intera*);
void Mesh_viewer_Arcroll_key_callback(Mesh_viewer_Intera*);
void Mesh_viewer_Arcroll_scroll_callback(Mesh_viewer_Intera*,double,double);
#ifdef __cplusplus
}
#endif

#endif
