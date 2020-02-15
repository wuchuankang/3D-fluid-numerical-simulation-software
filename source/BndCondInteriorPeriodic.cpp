

#include "Mesh.h"
#include "Geometry.h"

using namespace std;

void CMesh::BndCondInterior(CBlocks &A, CBlocks &B, const TInterior &Interior)								   
{
	int NDummy = CGeometry::NDummy;
	int nd, fc_a, fc_b, orient;
	int is_a, ie_a, js_a, je_a, ks_a, ke_a,
		is_b, ie_b, js_b, je_b, ks_b, ke_b;

	orient = Interior.orient;
	fc_a = Interior.Face_type_1;
	fc_b = Interior.Face_type_2;
	is_a = Interior.is_1;		is_b = Interior.is_2;
	ie_a = Interior.ie_1;		ie_b = Interior.ie_2;
	js_a = Interior.js_1;		js_b = Interior.js_2;
	je_a = Interior.je_1;		je_b = Interior.je_2;
	ks_a = Interior.ks_1;		ks_b = Interior.ks_2;
	ke_a = Interior.ke_1;		ke_b = Interior.ke_2;

	for(nd=1;nd<=NDummy;nd++)
	{
		ShadowOut(A, fc_a, nd, is_a, ie_a, js_a, je_a, ks_a, ke_a);
		Trans(orient,0.0);
		ShadowIn (B, fc_b, nd, is_b, ie_b, js_b, je_b, ks_b, ke_b);	
		ShadowOut(B, fc_b, nd, is_b, ie_b, js_b, je_b, ks_b, ke_b);
		Trans(orient,0.0);
		ShadowIn(A, fc_a, nd, is_a, ie_a, js_a, je_a, ks_a, ke_a);	
	}
}

void CMesh::BndCondPeriodic(CBlocks &A, CBlocks &B, const TPeriodic &Periodic)									
{
	int NDummy = CGeometry::NDummy;
	int nd, fc_a, fc_b, orient;
	int is_a, ie_a, js_a, je_a, ks_a, ke_a,
		is_b, ie_b, js_b, je_b, ks_b, ke_b;
	REAL th_1, th_2;

	orient = Periodic.orient;
	fc_a = Periodic.Face_type_1;
	fc_b = Periodic.Face_type_2;
	is_a = Periodic.is_1;		is_b = Periodic.is_2;
	ie_a = Periodic.ie_1;		ie_b = Periodic.ie_2;
	js_a = Periodic.js_1;		js_b = Periodic.js_2;
	je_a = Periodic.je_1;		je_b = Periodic.je_2;
	ks_a = Periodic.ks_1;		ks_b = Periodic.ks_2;
	ke_a = Periodic.ke_1;		ke_b = Periodic.ke_2;

	th_1 = ATAN(A.Geometry.Node_coord[is_a][js_a][ks_a].y / A.Geometry.Node_coord[is_a][js_a][ks_a].x);
	th_2 = ATAN(B.Geometry.Node_coord[is_b][js_b][ks_b].y / B.Geometry.Node_coord[is_b][js_b][ks_b].x);
	if((th_2 - th_1)<0.0)
	{
		A.Geometry.delta_th = 0.0f-A.Geometry.delta_th;
	}
	for(nd=1;nd<=NDummy;nd++)
	{
		ShadowOut(A, fc_a, nd, is_a, ie_a, js_a, je_a, ks_a, ke_a);
		Trans(orient,A.Geometry.delta_th);
		ShadowIn (B, fc_b, nd, is_b, ie_b, js_b, je_b, ks_b, ke_b);
		ShadowOut(B, fc_b, nd, is_b, ie_b, js_b, je_b, ks_b, ke_b);
		Trans(orient, 0.0f-A.Geometry.delta_th);
		ShadowIn(A, fc_a, nd, is_a, ie_a, js_a, je_a, ks_a, ke_a);
	}
}