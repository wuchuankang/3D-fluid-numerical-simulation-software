

#include "defs.h"

// constructor function
TNode::TNode() { x=0.0; y=0.0; z=0.0; }

TFace::TFace() 
{ 
	r=0.0;          th=0.0;			s=0.0;
	gradx_kesi=0.0; grady_kesi=0.0; gradz_kesi=0.0;
	gradx_yita=0.0; grady_yita=0.0;	gradz_yita=0.0;
	gradx_zita=0.0; grady_zita=0.0;	gradz_zita=0.0;
}

TCell::TCell()
{
	dis_Wall=0.0;
	r=0.0;			th=0.0;			Vol=0.0;	
	gradx_kesi=0.0;	grady_kesi=0.0; gradz_kesi=0.0;
	gradx_yita=0.0;	grady_yita=0.0; gradz_yita=0.0;
	gradx_zita=0.0;	grady_zita=0.0; gradz_zita=0.0;	
}


TConsVars::TConsVars() { dens=0.0; xmom=0.0; ymom=0.0; zmom=0.0; ener=0.0; }
	
TPrimVars::TPrimVars() { dens=0.0; uvel=0.0; vvel=0.0; wvel=0.0; press=0.0; temp=0.0; }

TVisVars::TVisVars()   { cond=0.0; mue=0.0; }

TEddyVars::TEddyVars() { ev=0.0; }

// override operators
// about TConsVars
TConsVars operator+ (const TConsVars &A, const TConsVars &B)
{
	TConsVars C;
	C.dens = A.dens + B.dens;
	C.xmom = A.xmom + B.xmom;
	C.ymom = A.ymom + B.ymom;
	C.zmom = A.zmom + B.zmom;
	C.ener = A.ener + B.ener;
	return C;
}

TConsVars operator+ (const TConsVars &A, REAL b)
{
	TConsVars C;
	C.dens = A.dens + b;
	C.xmom = A.xmom + b;
	C.ymom = A.ymom + b;
	C.zmom = A.zmom + b;
	C.ener = A.ener + b;
	return C;
}
TConsVars operator+ (REAL b, const TConsVars &A)
{
	TConsVars C;
	C.dens = A.dens + b;
	C.xmom = A.xmom + b;
	C.ymom = A.ymom + b;
	C.zmom = A.zmom + b;
	C.ener = A.ener + b;
	return C;
}

TConsVars operator- (const TConsVars &A, const TConsVars &B)
{
	TConsVars C;
	C.dens = A.dens - B.dens;
	C.xmom = A.xmom - B.xmom;
	C.ymom = A.ymom - B.ymom;
	C.zmom = A.zmom - B.zmom;
	C.ener = A.ener - B.ener;
	return C;
}

TConsVars operator- (const TConsVars &A, REAL b)
{
	TConsVars C;
	C.dens = A.dens - b;
	C.xmom = A.xmom - b;
	C.ymom = A.ymom - b;
	C.zmom = A.zmom - b;
	C.ener = A.ener - b;
	return C;
}
TConsVars operator- (REAL b, const TConsVars &A)
{
	TConsVars C;
	C.dens = b-A.dens;
	C.xmom = b-A.xmom;
	C.ymom = b-A.ymom;
	C.zmom = b-A.zmom;
	C.ener = b-A.ener;
	return C;
}
TConsVars operator* (const TConsVars &A, const TConsVars &B)
{
	TConsVars C;
	C.dens = A.dens*B.dens;
	C.xmom = A.xmom*B.xmom;
	C.ymom = A.ymom*B.ymom;
	C.zmom = A.zmom*B.zmom;
	C.ener = A.ener*B.ener;
	return C;
}
TConsVars operator* (const TConsVars &A, REAL b)
{
	TConsVars C;
	C.dens = A.dens*b;
	C.xmom = A.xmom*b;
	C.ymom = A.ymom*b;
	C.zmom = A.zmom*b;
	C.ener = A.ener*b;
	return C;
}
TConsVars operator* (REAL b, const TConsVars &A)
{
	TConsVars C;
	C.dens = A.dens*b;
	C.xmom = A.xmom*b;
	C.ymom = A.ymom*b;
	C.zmom = A.zmom*b;
	C.ener = A.ener*b;
	return C;
}	

TConsVars operator/ (const TConsVars &A, const TConsVars &B)
{
	TConsVars C;
	C.dens = A.dens/B.dens;
	C.xmom = A.xmom/B.xmom;
	C.ymom = A.ymom/B.ymom;
	C.zmom = A.zmom/B.zmom;
	C.ener = A.ener/B.ener;
	return C;
}

TConsVars operator/ (const TConsVars &A, REAL b )
{
	TConsVars C;
	C.dens = A.dens/b;
	C.xmom = A.xmom/b;
	C.ymom = A.ymom/b;
	C.zmom = A.zmom/b;
	C.ener = A.ener/b;
	return C;
}

TConsVars operator/ (REAL b,  const TConsVars &A)
{
	TConsVars C;
	C.dens = b/A.dens;
	C.xmom = b/A.xmom;
	C.ymom = b/A.ymom;
	C.zmom = b/A.zmom;
	C.ener = b/A.ener;
	return C;
}

TConsVars POW(const TConsVars &A, REAL b )
{
	TConsVars C;
	C.dens = POW(A.dens, b);
	C.xmom = POW(A.xmom, b);
	C.ymom = POW(A.ymom, b);
	C.zmom = POW(A.zmom, b);
	C.ener = POW(A.ener, b);
	return C;
}

// about TPrimVars
TPrimVars operator+ (const TPrimVars &A, const TPrimVars &B)
{
	TPrimVars C;
	C.dens  = A.dens  + B.dens;
	C.uvel  = A.uvel  + B.uvel;
	C.vvel  = A.vvel  + B.vvel;
	C.wvel  = A.wvel  + B.wvel;
	C.temp  = A.temp  + B.temp;
	C.press = A.press + B.press;
	return C;
}

TPrimVars operator+ (const TPrimVars &A, REAL b)
{
	TPrimVars C;
	C.dens  = A.dens  + b;
	C.uvel  = A.uvel  + b;
	C.vvel  = A.vvel  + b;
	C.wvel  = A.wvel  + b;
	C.temp  = A.temp  + b;
	C.press = A.press + b;
	return C;
}

TPrimVars operator+ (REAL b, const TPrimVars &A)
{
	TPrimVars C;
	C.dens  = A.dens  + b;
	C.uvel  = A.uvel  + b;
	C.vvel  = A.vvel  + b;
	C.wvel  = A.wvel  + b;
	C.temp  = A.temp  + b;
	C.press = A.press + b;
	return C;
}

TPrimVars operator- (const TPrimVars &A, const TPrimVars &B)
{
	TPrimVars C;
	C.dens  = A.dens  - B.dens;
	C.uvel  = A.uvel  - B.uvel;
	C.vvel  = A.vvel  - B.vvel;
	C.wvel  = A.wvel  - B.wvel;
	C.temp  = A.temp  - B.temp;
	C.press = A.press - B.press;
	return C;
}

TPrimVars operator- (REAL b, const TPrimVars &A)
{
	TPrimVars C;
	C.dens  = b - A.dens;
	C.uvel  = b - A.uvel;
	C.vvel  = b - A.vvel;
	C.wvel  = b - A.wvel;
	C.temp  = b - A.temp;
	C.press = b - A.press;
	return C;
}

TPrimVars operator* (const TPrimVars &A, const TPrimVars &B)
{
	TPrimVars C;
	C.dens  = A.dens * B.dens;
	C.uvel  = A.uvel * B.uvel;
	C.vvel  = A.vvel * B.vvel;
	C.wvel  = A.wvel * B.wvel;
	C.temp  = A.temp * B.temp;
	C.press = A.press * B.press;
	return C;
}

TPrimVars operator* (REAL b, const TPrimVars &A)
{
	TPrimVars C;
	C.dens  = b * A.dens;
	C.uvel  = b * A.uvel;
	C.vvel  = b * A.vvel;
	C.wvel  = b * A.wvel;
	C.temp  = b * A.temp;
	C.press = b * A.press;
	return C;
}

TPrimVars operator/ (const TPrimVars &A, REAL b)
{
	TPrimVars C;
	C.dens  = A.dens / b;
	C.uvel  = A.uvel / b;
	C.vvel  = A.vvel / b;
	C.wvel  = A.wvel / b;
	C.temp  = A.temp / b;
	C.press = A.press/ b;
	return C;
}

TPrimVars operator/ (REAL b, const TPrimVars &A)
{
	TPrimVars C;
	C.dens  = b / A.dens ;
	C.uvel  = b / A.uvel ;
	C.vvel  = b / A.vvel ;
	C.wvel  = b / A.wvel ;
	C.temp  = b / A.temp ;
	C.press = b / A.press;
	return C;
}

TPrimVars operator/ (const TPrimVars &A, const TPrimVars &B)
{
	TPrimVars C;
	C.dens  = A.dens / B.dens ;
	C.uvel  = A.uvel / B.uvel ;
	C.vvel  = A.vvel / B.vvel ;
	C.wvel  = A.wvel / B.wvel ;
	C.temp  = A.temp / B.temp ;
	C.press = A.press/ B.press;
	return C;
}

TPrimVars POW(const TPrimVars &A, REAL b )
{
	TPrimVars C;
	C.dens  = POW(A.dens, b);
	C.uvel  = POW(A.uvel, b);
	C.vvel  = POW(A.vvel, b);
	C.wvel  = POW(A.wvel, b);
	C.temp  = POW(A.temp, b);
	C.press = POW(A.press, b);
	return C;
}

// about TNode
REAL operator* (const TNode &A, const TNode &B)
{
	REAL C = A.x*B.x + A.y*B.y + A.z*B.z;
	return C;
}

TNode operator- (const TNode &A, const TNode &B)
{
	TNode C;
	C.x = A.x-B.x;
	C.y = A.y-B.y;
	C.z = A.z-B.z;
	return C;
}

TNode operator+ (const TNode &A, const TNode &B)
{
	TNode C;
	C.x = A.x+B.x;
	C.y = A.y+B.y;
	C.z = A.z+B.z;
	return C;
}

TNode operator* (const TNode &A, REAL b)
{
	TNode C;
	C.x = A.x*b;
	C.y = A.y*b;
	C.z = A.z*b;
	return C;
}

TNode operator* (REAL b, const TNode &A)
{
	TNode C;
	C.x = A.x*b;
	C.y = A.y*b;
	C.z = A.z*b;
	return C;
}

TNode operator/ (const TNode &A, REAL b)
{
	TNode C;
	C.x = A.x/b;
	C.y = A.y/b;
	C.z = A.z/b;
	return C;
}

TConsVars operator* (REAL A[][6], TConsVars &b)
{	
	TConsVars c; 

	c.dens = A[1][1]*b.dens + A[1][2]*b.xmom + A[1][3]*b.ymom + A[1][4]*b.zmom + A[1][5]*b.ener;
   	c.xmom = A[2][1]*b.dens + A[2][2]*b.xmom + A[2][3]*b.ymom + A[2][4]*b.zmom + A[2][5]*b.ener;
	c.ymom = A[3][1]*b.dens + A[3][2]*b.xmom + A[3][3]*b.ymom + A[3][4]*b.zmom + A[3][5]*b.ener;
	c.zmom = A[4][1]*b.dens + A[4][2]*b.xmom + A[4][3]*b.ymom + A[4][4]*b.zmom + A[4][5]*b.ener;
	c.ener = A[5][1]*b.dens + A[5][2]*b.xmom + A[5][3]*b.ymom + A[5][4]*b.zmom + A[5][5]*b.ener;

	return c;
}


REAL MAX4(REAL a1, REAL a2, REAL a3, REAL a4)
{
	REAL max;
	max=a1;
	if(max<a2)	max=a2;
	if(max<a3)	max=a3;
	if(max<a4)	max=a4;	
	return (max);
}

REAL MAX3(REAL a1, REAL a2, REAL a3)
{
	REAL max;
	max=a1;
	if(max<a2)	max=a2;
	if(max<a3)	max=a3;
	return (max);
}

int MAX3(int a1, int a2, int a3)
{
	int max;
	max=a1;
	if(max<a2)	max=a2;
	if(max<a3)	max=a3;
	return (max);
}

REAL MIN3(REAL a1, REAL a2, REAL a3)
{
	REAL min;
	min=a1;
	if(min>a2)	min=a2;
	if(min>a3)	min=a3;
	return (min);
}

int MIN3(int a1, int a2, int a3)
{
	int min;
	min=a1;
	if(min>a2)	min=a2;
	if(min>a3)	min=a3;
	return (min);
}

int NT;