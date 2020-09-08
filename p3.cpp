#include		<Windows.h>
//#define			min(x,y) (x<y?x:y)
//#define			max(x,y) (x>y?x:y)
#include		<sys\stat.h>
#include		<fstream>
#include		<string>
#include		<vector>
#include		<algorithm>
#include		<functional>
#include		<Shlwapi.h>
#include		<gdiplus.h>

#include		<gl/gl.h>//
#include		<gl/glu.h>//
//#include		<glut.h>
//#include		<glext.h>
#include		"openglut.h"

//#define			QUICKHULL_IMPLEMENTATION
//#include		"quickhull.h"
#pragma			comment(lib, "Shlwapi.lib")
#pragma			comment(lib, "gdiplus.lib")
#pragma			comment(lib, "OpenGL32.lib")
#pragma			comment(lib, "GLU32.lib")
//#pragma			comment(lib, "GLUT32.lib")
#pragma			comment(lib, "OpenGLUT.lib")
void			(__stdcall *glWindowPos2i)(int, int)=nullptr;
//PFNGLWINDOWPOS2IPROC glWindowPos2i;

//	#define		AA_MOUSEMOVE_TRIGGER//draw on any mouse move for aa

#define			MOUSE_SENSITIVITY	.003f
#define			DTANFOV				1.1f
#define			DFOV				1.1f
#define			TIMER_ELAPSE		10
int				nRandPols=5, nRandTPols=5;//5	100

unsigned		EPA_iteration=0, n_EPA_iterations=0, n_GJK_iterations=0;
int				scenario_number=0, frame_number=0,
				n_iterations=10;

float			timescale=0.5,//0.05f
				gravity=0.98f,//9.8f	0.098f
				zground=0,
				xleft=0, yfront=0,
				beta=0.2f,//bias factor [0.1 0.3]
				C_R=0.3f,//coefficient of restitution
				boundary_friction=0.5f,
				object_friction=0.5f,
				proximity_limit=5,
				penetration_slop=0.5;

bool			pause=true, tick=true, info=true,
				accumulated_impulses=true,
			//	warm_starting=true,
				warm_starting=false,
				lsw_transparency_multiply=false;
const float		_pi=acos(-1.f), _2pi=2*_pi, pi_2=_pi*0.5f, sqrt2=sqrt(2.f);
#if __cplusplus<201103
#pragma			warning(push)
#pragma			warning(disable:4723)
const float		zero=0, infinity=1/zero;
#pragma			warning(pop)
inline bool		isinf(float x){return ((int&)x&0x7FFFFFFF)==(int&)infinity;}
inline bool		isnan(float x){return x!=x;}
#endif
inline bool		is_nan_or_inf(float x){return (int&)x==0x7F800000;}
inline __m128	_mm_sin_ps(__m128 angle)
{
	const __m128 m_pi=_mm_set1_ps(_pi), m_2pi=_mm_set1_ps(_2pi);
	const __m128 half=_mm_set1_ps(0.5f);
	const __m128 f3=_mm_set1_ps(6), f5=_mm_set1_ps(120), f7=_mm_set1_ps(5040), f9=_mm_set1_ps(362880), f11=_mm_set1_ps(39916800);

	//bring angle to the range [-pi, pi[
	//angle -= floor(angle/2pi)*2pi		X [0, 2pi[
	//angle -= floor(angle/2pi-0.5)*2pi - 2pi		[-pi, pi[
	__m128 angle_pi=_mm_div_ps(angle, m_pi);
	angle_pi=_mm_sub_ps(angle_pi, half);
	angle_pi=_mm_floor_ps(angle_pi);
	angle_pi=_mm_mul_ps(angle_pi, m_2pi);
	angle=_mm_sub_ps(angle, angle_pi);
	angle=_mm_sub_ps(angle, m_2pi);

	__m128 m_sinx=angle;
	__m128 angle2=_mm_mul_ps(angle, angle);
	__m128 angle_n=_mm_mul_ps(angle, angle2);//x^3
	__m128 term=_mm_div_ps(angle_n, f3);
	m_sinx=_mm_sub_ps(m_sinx, term);

	angle_n=_mm_mul_ps(angle_n, angle2);//x^5
	term=_mm_div_ps(angle_n, f5);
	m_sinx=_mm_add_ps(m_sinx, term);

	angle_n=_mm_mul_ps(angle_n, angle2);//x^7
	term=_mm_div_ps(angle_n, f7);
	m_sinx=_mm_sub_ps(m_sinx, term);

	angle_n=_mm_mul_ps(angle_n, angle2);//x^9
	term=_mm_div_ps(angle_n, f9);
	m_sinx=_mm_add_ps(m_sinx, term);

	angle_n=_mm_mul_ps(angle_n, angle2);//x^11
	term=_mm_div_ps(angle_n, f11);
	m_sinx=_mm_sub_ps(m_sinx, term);
	return m_sinx;
}
inline __m128	_mm_cos_ps(__m128 angle)
{
	const __m128 m_pi_2=_mm_set1_ps(_pi/2);
	return _mm_sin_ps(_mm_sub_ps(m_pi_2, angle));
}
void			render				();
struct			Vector2f
{
	float x, y;
	Vector2f():x(0), y(0){}
	Vector2f(float x, float y):x(x), y(y){}
	void set(float x, float y){this->x=x, this->y=y;}
	Vector2f& operator+=(Vector2f const &b){x+=b.x, y+=b.y; return *this;}
	Vector2f& operator-=(Vector2f const &b){x-=b.x, y-=b.y; return *this;}
	Vector2f& operator+=(float x){this->x+=x, y+=x; return *this;}
	Vector2f& operator-=(float x){this->x-=x, y-=x; return *this;}
	Vector2f& operator*=(float x){this->x*=x, y*=x; return *this;}
	Vector2f& operator/=(float x){this->x/=x, y/=x; return *this;}
	float dot(Vector2f const &other)const{return x*other.x+y*other.y;}
	float cross(Vector2f const &other)const{return x*other.y-y*other.x;}
	float magnitude()const{return sqrt(x*x+y*y);}
	float mag_sq()const{return x*x+y*y;}
	float angle()const{return atan(y/x);}
	float angle2()const{return atan2(y, x);}
};
inline Vector2f	operator*(Vector2f const &p, float x){return Vector2f(p.x*x, p.y*x);}
inline Vector2f	operator*(float x, Vector2f const &p){return Vector2f(p.x*x, p.y*x);}
inline Vector2f	operator/(Vector2f const &p, float x){return Vector2f(p.x/x, p.y/x);}
inline Vector2f	operator+(Vector2f const &a, Vector2f const &b){return Vector2f(a.x+b.x, a.y+b.y);}
inline Vector2f	operator-(Vector2f const &a, Vector2f const &b){return Vector2f(a.x-b.x, a.y-b.y);}
inline Vector2f	operator+(Vector2f const &p, float x){return Vector2f(p.x+x, p.y+x);}
inline Vector2f	operator-(Vector2f const &p, float x){return Vector2f(p.x-x, p.y-x);}
inline bool		operator==(Vector2f const &a, Vector2f const &b){return a.x==b.x&&a.y==b.y;}
inline bool		operator!=(Vector2f const &a, Vector2f const &b){return a.x!=b.x||a.y!=b.y;}
inline Vector2f	operator-(Vector2f const &p){return Vector2f(-p.x, -p.y);}
struct			Matrix22
{
	float a, b,//row 0
		c, d;//row 1
	Matrix22():a(0), b(0), c(0), d(0){}
	Matrix22(float a, float b, float c, float d):a(a), b(b), c(c), d(d){}
};
Vector2f		operator*(Matrix22 const &m, Vector2f const &v){return Vector2f(m.a*v.x+m.b*v.y, m.c*v.x+m.d*v.y);}
struct			Vector3f
{
	float x, y, z;
	Vector3f():x(0), y(0), z(0){}
	Vector3f(float x, float y, float z):x(x), y(y), z(z){}
	void set(float x, float y, float z){this->x=x, this->y=y, this->z=z;}
	Vector3f& operator+=(Vector3f const &b){x+=b.x, y+=b.y, z+=b.z; return *this;}
	Vector3f& operator-=(Vector3f const &b){x-=b.x, y-=b.y, z-=b.z; return *this;}
	Vector3f& operator+=(float x){this->x+=x, y+=x, z+=x; return *this;}
	Vector3f& operator-=(float x){this->x-=x, y-=x, z-=x; return *this;}
	Vector3f& operator*=(float x){this->x*=x, y*=x, z*=x; return *this;}
	Vector3f& operator/=(float x){this->x/=x, y/=x, z/=x; return *this;}
	float dot(Vector3f const &other)const{return x*other.x+y*other.y+z*other.z;}
	Vector3f cross(Vector3f const &other)const{return Vector3f(y*other.z-z*other.y, z*other.x-x*other.z, x*other.y-y*other.x);}
	Vector3f triple_product(Vector3f const &b, Vector3f const &c)const;
	float magnitude()const{return sqrt(x*x+y*y+z*z);}
	float mag_sq()const{return x*x+y*y+z*z;}
	float theta()const{return atan(z/sqrt(x*x+y*y));}//vertical angle
	//float theta2()const{return atan2(y, x);}
	float phi(){return atan(y/x);}//horizontal angle
	float phi2(){return atan2(y, x);}
	bool isnan(){return x!=x||y!=y||z!=z;}
	bool isnan_or_inf(){return x!=x||y!=y||z!=z||abs(x)==infinity||abs(y)==infinity||abs(z)==infinity;}
};
inline Vector3f	operator*(Vector3f const &p, float x){return Vector3f(p.x*x, p.y*x, p.z*x);}
inline Vector3f	operator*(float x, Vector3f const &p){return Vector3f(p.x*x, p.y*x, p.z*x);}
inline Vector3f	operator/(Vector3f const &p, float x){return Vector3f(p.x/x, p.y/x, p.z/x);}
inline Vector3f	operator+(Vector3f const &a, Vector3f const &b){return Vector3f(a.x+b.x, a.y+b.y, a.z+b.z);}
inline Vector3f	operator-(Vector3f const &a, Vector3f const &b){return Vector3f(a.x-b.x, a.y-b.y, a.z-b.z);}
inline Vector3f	operator+(Vector3f const &p, float x){return Vector3f(p.x+x, p.y+x, p.z+x);}
inline Vector3f	operator-(Vector3f const &p, float x){return Vector3f(p.x-x, p.y-x, p.z-x);}
inline bool		operator==(Vector3f const &a, Vector3f const &b){return a.x==b.x&&a.y==b.y&&a.z==b.z;}
inline bool		operator!=(Vector3f const &a, Vector3f const &b){return a.x!=b.x||a.y!=b.y||a.z!=b.z;}
inline Vector3f	operator-(Vector3f const &p){return Vector3f(-p.x, -p.y, -p.z);}
Vector3f		Vector3f::triple_product(Vector3f const &b, Vector3f const &c)const{return this->dot(c)*b-this->dot(b)*c;}
struct			Triangle0
{
	Vector3f p1, p2, p3;
	Vector2f tx1, tx2, tx3;
	Matrix22 txm;
//	float x1, y1, z1, x2, y2, z2, x3, y3, z3, tx1x, tx1y, tx2x, tx2y, tx3x, tx3y, txXX, txXY, txYX, txYY;
	char tx_idx;
	int color;
} *env=0, *tenv=0;
float			adjust_angle(float angle)//reduce angle to the range [0, 2pi[
{
	float angle_pi2=angle/_2pi;
	angle=_2pi*(angle_pi2-floor(angle_pi2));
	if(angle<0||angle>=_2pi)
		return 0;
	return angle;
}
//struct			HalfEdge
//{
//	int prev, next;
//};
struct			Triangle
{
	int a, b, c;//point indeces
	bool textured;
	int tc;//texture index / color
	Vector2f tx1, tx2, tx3;//coordinates of vertices on texture
	Matrix22 txm;
	Triangle():a(-1), b(-1), c(-1), textured(false){}
	Triangle(int a, int b, int c, int color):a(a), b(b), c(c), textured(false), tc(color){}
	Triangle(int a, int b, int c, int tx_idx, Vector2f const &tx1, Vector2f const &tx2, Vector2f const &tx3):a(a), b(b), c(c), textured(true), tx1(tx1), tx2(tx2), tx3(tx3)
	{
		//set texture matrix
	}
	void set(int a, int b, int c, int tx_idx, Vector2f const &tx1, Vector2f const &tx2, Vector2f const &tx3)
	{
		this->a=a, this->b=b, this->c=c, textured=true, this->tx1=tx1, this->tx2=tx2, this->tx3=tx3;
		//set texture matrix
	}
	void set(int a, int b, int c, int color){this->a=a, this->b=b, this->c=c, textured=false, tc=color;}
};
struct			Edge
{
	int v1, v2, tr1, tr2;
	Edge():v1(-1), v2(-1), tr1(-1), tr2(-1){}
	Edge(int v1, int v2, int tr1, int tr2):v1(v1), v2(v2), tr1(tr1), tr2(tr2){}
	void set(int v1, int v2, int tr1, int tr2){this->v1=v1, this->v2=v2, this->tr1=tr1, this->tr2=tr2;}
};
struct			Matrix33
{
	float
		a11, a12, a13,
		a21, a22, a23,
		a31, a32, a33;
	Matrix33():a11(0), a12(0), a13(0), a21(0), a22(0), a23(0), a31(0), a32(0), a33(0){}
	Matrix33(float a11, float a12, float a13, float a21, float a22, float a23, float a31, float a32, float a33):a11(a11), a12(a12), a13(a13), a21(a21), a22(a22), a23(a23), a31(a31), a32(a32), a33(a33){}
	void transpose_this()
	{
		std::swap(a21, a12);
		std::swap(a31, a13);
		std::swap(a32, a23);
	}
	Matrix33 transpose()
	{
		Matrix33 m(
			a11, a21, a31,
			a12, a22, a32,
			a13, a23, a33);
		return m;
	}
	float determinant()
	{
		return a11*(a22*a33-a23*a32)-a12*(a21*a33-a23*a31)+a13*(a21*a32-a22*a31);
	}
	float determinant_diagonal(){return a11*a22*a33;}
	Matrix33 inverse()
	{
		float det=determinant();
		det=1/det;
		Matrix33 m(
			det*(a22*a33-a23*a32), det*(a13*a32-a12*a33), det*(a12*a23-a13*a22),
			det*(a23*a31-a21*a33), det*(a11*a33-a13*a31), det*(a13*a21-a11*a23),
			det*(a21*a32-a22*a31), det*(a12*a31-a11*a32), det*(a11*a22-a12*a21));
		return m;
	}
	Matrix33 inverse_diagonal()
	{
		return Matrix33(
			1/a11, 0, 0,
			0, 1/a22, 0,
			0, 0, 1/a33);
	}
	Vector3f multiply_diagonal(Vector3f const &v){return Vector3f(a11*v.x, a22*v.y, a33*v.z);}
	//bool inverse(Matrix33 &m)
	//{
	//	float det=a11*(a22*a33-a23*a32)-a12*(a21*a33-a23*a31)+a13*(a21*a32-a22*a31);
	//	if(det)
	//	{
	//		det=1/det;
	//		m.a11=det*(a22*a33-a23*a32), m.a12=det*(a13*a32-a12*a33), m.a13=det*(a12*a23-a13*a22);
	//		m.a21=det*(a23*a31-a21*a33), m.a22=det*(a11*a33-a13*a31), m.a23=det*(a13*a21-a11*a23);
	//		m.a31=det*(a21*a32-a22*a31), m.a32=det*(a12*a31-a11*a32), m.a33=det*(a11*a22-a12*a21);
	//		return true;
	//	}
	//	return false;
	//}
};
Vector3f		operator*(Matrix33 const &m, Vector3f const &v)
{
	return Vector3f(
		m.a11*v.x+m.a12*v.y+m.a13*v.z,
		m.a21*v.x+m.a22*v.y+m.a23*v.z,
		m.a31*v.x+m.a32*v.y+m.a33*v.z);
}
struct			Quaternion
{
	float w, x, y, z;
	Quaternion():w(0), x(0), y(0), z(0){}
	Quaternion(float w, float x, float y, float z):w(w), x(x), y(y), z(z){}
	void set(float w, float x, float y, float z){this->w=w, this->x=x, this->y=y, this->z=z;}
	Quaternion conjugate(){return Quaternion(w, -x, -y, -z);}
	float magnitude(){return sqrt(w*w+x*x+y*y+z*z);}
	Vector3f mul_pure(Quaternion const &q2)
	{
		const __m128 mask_pnnn=_mm_castsi128_ps(_mm_set_epi32(0x80000000, 0x80000000, 0x80000000, 0)),
			mask_pppn=_mm_castsi128_ps(_mm_set_epi32(0x80000000, 0, 0, 0)),
			mask_pnpp=_mm_castsi128_ps(_mm_set_epi32(0, 0, 0x80000000, 0)),
			mask_ppnp=_mm_castsi128_ps(_mm_set_epi32(0, 0x80000000, 0, 0));
		__m128 wxyz1=_mm_loadu_ps(&w), wxyz2=_mm_loadu_ps(&q2.w);
		__m128 xwzy2=_mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(wxyz2), _MM_SHUFFLE(2, 3, 0, 1)));
		__m128 r_x=_mm_mul_ps(wxyz1, xwzy2);
		__m128 yzwx2=_mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(wxyz2), _MM_SHUFFLE(1, 0, 3, 2)));
		__m128 r_y=_mm_mul_ps(wxyz1, yzwx2);
		__m128 zyxw2=_mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(wxyz2), _MM_SHUFFLE(0, 1, 2, 3)));
		__m128 r_z=_mm_mul_ps(wxyz1, zyxw2);
		r_x=_mm_xor_ps(r_x, mask_pppn);
		r_y=_mm_xor_ps(r_y, mask_pnpp);
		r_z=_mm_xor_ps(r_z, mask_ppnp);
		__m128 r_wx=_mm_hadd_ps(_mm_setzero_ps(), r_x);
		__m128 r_yz=_mm_hadd_ps(r_y, r_z);
		__m128 r_wxyz=_mm_hadd_ps(r_wx, r_yz);
		Vector3f r;
		r.x=r_wxyz.m128_f32[1], r.y=r_wxyz.m128_f32[2], r.z=r_wxyz.m128_f32[3];
		return r;
		//return Vector3f(
		//	w*q2.x+x*q2.w+y*q2.z-z*q2.y,
		//	w*q2.y-x*q2.z+y*q2.w+z*q2.x,
		//	w*q2.z+x*q2.y-y*q2.x+z*q2.w);
	}
};
Quaternion	 operator*(Quaternion const &q1, Quaternion const &q2)//q2 followed by q1
{
	const __m128 mask_pnnn=_mm_castsi128_ps(_mm_set_epi32(0x80000000, 0x80000000, 0x80000000, 0)),
		mask_pppn=_mm_castsi128_ps(_mm_set_epi32(0x80000000, 0, 0, 0)),
		mask_pnpp=_mm_castsi128_ps(_mm_set_epi32(0, 0, 0x80000000, 0)),
		mask_ppnp=_mm_castsi128_ps(_mm_set_epi32(0, 0x80000000, 0, 0));
	__m128 wxyz1=_mm_loadu_ps(&q1.w), wxyz2=_mm_loadu_ps(&q2.w);
	__m128 r_w=_mm_mul_ps(wxyz1, wxyz2);
	__m128 xwzy2=_mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(wxyz2), _MM_SHUFFLE(2, 3, 0, 1)));//2, 3, 0, 1
	__m128 r_x=_mm_mul_ps(wxyz1, xwzy2);
	__m128 yzwx2=_mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(wxyz2), _MM_SHUFFLE(1, 0, 3, 2)));//1, 0, 3, 2
	__m128 r_y=_mm_mul_ps(wxyz1, yzwx2);
	__m128 zyxw2=_mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(wxyz2), _MM_SHUFFLE(0, 1, 2, 3)));//0, 1, 2, 3
	__m128 r_z=_mm_mul_ps(wxyz1, zyxw2);
	r_w=_mm_xor_ps(r_w, mask_pnnn);
	r_x=_mm_xor_ps(r_x, mask_pppn);
	r_y=_mm_xor_ps(r_y, mask_pnpp);
	r_z=_mm_xor_ps(r_z, mask_ppnp);
	__m128 r_wx=_mm_hadd_ps(r_w, r_x);
	__m128 r_yz=_mm_hadd_ps(r_y, r_z);
	__m128 r_wxyz=_mm_hadd_ps(r_wx, r_yz);
	Quaternion r;
	_mm_storeu_ps(&r.w, r_wxyz);
	return r;
	//return Quaternion(
	//	q1.w*q2.w-q1.x*q2.x-q1.y*q2.y-q1.z*q2.z,
	//	q1.w*q2.x+q1.x*q2.w+q1.y*q2.z-q1.z*q2.y,
	//	q1.w*q2.y-q1.x*q2.z+q1.y*q2.w+q1.z*q2.x,
	//	q1.w*q2.z+q1.x*q2.y-q1.y*q2.x+q1.z*q2.w);
}
struct			PhysObject
{
	int idx;
	std::vector<Vector3f> w,//points in world coordinates
		vv;//points in local coordinates
	std::vector<Triangle> tr;
	std::vector<Edge> e;
	float r;//reach (distance to farthest point from centroid)

	float density,//properties
		mass, inv_mass;
	Matrix33 I, inv_I;
	bool I_diagonal;
	float friction;

	Vector3f p;//position (centroid)
	float yaw, pitch, roll;//orientation: yaw: z, pitch: y, roll: x
	Matrix33 o, inv_o;//orientation
	Quaternion or;
	Vector3f vr;//wx, wy, wz
	//float wx, wy, wz;//angular velocities

	Vector3f v;//velocity vector
	void iterate(float delta)
	{
		if(frame_number==32)
			int LOL_1=0;
		p+=v*delta;

		float angle=vr.magnitude();
		Vector3f axis;
		if(angle<0.001f)//Taylor's expansion of sinc
			axis=vr*(0.5f*delta-delta*delta*delta*0.020833333333f*angle*angle);
		else
			axis=vr*sin(0.5f*angle*delta)/angle;
		Quaternion d_o(cos(0.5f*angle*delta), axis.x, axis.y, axis.z);
		or=d_o*or;

	//	Vector3f vrl=o*vr;
	//	yaw=adjust_angle(yaw+vrl.z*delta), pitch=adjust_angle(pitch+vrl.y*delta), roll=adjust_angle(roll+vrl.x*delta);

	//	yaw=adjust_angle(yaw+vr.z*delta), pitch=adjust_angle(pitch+vr.y*delta), roll=adjust_angle(roll+vr.x*delta);

		update_orientation(), update_world_coord();
	}
	void update_orientation()
	{
		float cy=cos(yaw), sy=sin(yaw), cp=cos(pitch), sp=sin(pitch), cr=cos(roll), sr=sin(roll);
		o.a11=cy*cp*cr-sy*sr, o.a12=-cy*cp*sr-sy*cr, o.a13=cy*sp;
		o.a21=sy*cp*cr+cy*sr, o.a22=-sy*cp*sr+cy*cr, o.a23=sy*sp;
		o.a31=-sp*cr, o.a32=sp*sr, o.a33=cp;
		inv_o=o.transpose();
	}
	void update_world_coord()
	{
		for(int k=0, kEnd=w.size();k<kEnd;++k)
		{
			Quaternion v(0, vv[k].x, vv[k].y, vv[k].z);
			w[k]=p+(or*v).mul_pure(or.conjugate());
			//v=or*v*or.conjugate();
			//w[k].set(v.x, v.y, v.z);
			//w[k]+=p;
		}
		//	w[k]=p+o*vv[k];
	}
	Vector3f world_to_local(Vector3f p)
	{
		p-=this->p;
		Quaternion v(0, p.x, p.y, p.z);
		return (or.conjugate()*v).mul_pure(or);
	//	return inv_o*p;
	}
	Vector3f local_to_world(Vector3f const &p)
	{
		Quaternion v(0, p.x, p.y, p.z);
		return this->p+(or*v).mul_pure(or.conjugate());
	//	return this->p+o*p;
	}
	Vector3f local_to_world_rotate_only(Vector3f const &p)
	{
		Quaternion v(0, p.x, p.y, p.z);
		return (or*v).mul_pure(or.conjugate());
	//	return o*p;
	}
	void set_position(float x, float y, float z, float yaw, float pitch, float roll)//thz, thy, thx
	{
		p.set(x, y, z);
		this->yaw=yaw, this->pitch=pitch, this->roll=roll;
		or.set(0, 1, 0, 0);
		or=Quaternion(cos(roll *0.5f), sin(roll *0.5f), 0, 0)*or;
		or=Quaternion(cos(pitch*0.5f), 0, sin(pitch*0.5f), 0)*or;
		or=Quaternion(cos(yaw  *0.5f), 0, 0, sin(yaw  *0.5f))*or;
	//	update_orientation();
	//	update_world_coord();
	}
	void set_velocity(float vx, float vy, float vz, float wx, float wy, float wz){v.set(vx, vy, vz), vr.x=wx, vr.y=wy, vr.z=wz;}

	Vector3f get_normal(int tr)const
	{
		auto &t=this->tr[tr];
		return (w[t.b]-w[t.a]).cross(w[t.c]-w[t.a]);//
	}

	//cuboid
	void set_properties_cuboid(float W, float L, float H, float density, int color)//dx, dy, dz
	{
		if(W<=0)
			W=1;
		if(L<=0)
			L=1;
		if(H<=0)
			H=1;
		vv.resize(8), w.resize(8), tr.resize(12);
		e.resize(18);
		//e.resize(12);
		float W2=W*0.5f, H2=H*0.5f, L2=L*0.5f;
		vv[0].set( W2,  L2,  H2);
		vv[1].set( W2,  L2, -H2);
		vv[2].set( W2, -L2, -H2);
		vv[3].set( W2, -L2,  H2);
		vv[4].set(-W2, -L2,  H2);
		vv[5].set(-W2, -L2, -H2);
		vv[6].set(-W2,  L2, -H2);
		vv[7].set(-W2,  L2,  H2);
		tr[0].set(0, 7, 3, color), tr[1].set(7, 4, 3, color);//top		winding: the normal (v2-v1).cross(v3-v1) points out of the object
		tr[2].set(1, 2, 6, color), tr[3].set(2, 5, 6, color);//bottom
		tr[4].set(2, 3, 4, color), tr[5].set(4, 5, 2, color);//left
		tr[6].set(0, 1, 7, color), tr[7].set(1, 6, 7, color);//right
		tr[8].set(1, 0, 3, color), tr[9].set(3, 2, 1, color);//front
		tr[10].set(7, 6, 5, color), tr[11].set(5, 4, 7, color);//back

		e[0].set(0, 7,  0, 6),		e[1].set(7, 4,  1, 11),		e[2].set(4, 3,  1, 4),		e[3].set(3, 0,  0, 8),		e[4].set(3, 7,  1, 0);
		e[5].set(0, 1,  6, 8),		e[6].set(1, 7,  6, 7),		e[7].set(6, 7,  7, 10),		e[8].set(7, 5,  10, 11);
		e[9].set(4, 5,  5, 11),		e[10].set(2, 4,  4, 5),		e[11].set(2, 3,  4, 9),		e[12].set(1, 3,  8, 9);
		e[13].set(1, 6,  2, 7),		e[14].set(6, 5,  7, 10),	e[15].set(5, 2,  3, 5),		e[16].set(1, 2,  2, 9),		e[17].set(2, 6,  2, 3);

		//tr[0].set(0, 7, 3, rand()<<16|rand()), tr[1].set(7, 4, 3, rand()<<16|rand());//top
		//tr[2].set(1, 2, 6, rand()<<16|rand()), tr[3].set(2, 5, 6, rand()<<16|rand());//bottom
		//tr[4].set(2, 3, 4, rand()<<16|rand()), tr[5].set(4, 5, 2, rand()<<16|rand());//left
		//tr[6].set(0, 1, 7, rand()<<16|rand()), tr[7].set(1, 6, 7, rand()<<16|rand());//right
		//tr[8].set(1, 0, 3, rand()<<16|rand()), tr[9].set(3, 2, 1, rand()<<16|rand());//front
		//tr[10].set(7, 6, 5, rand()<<16|rand()), tr[11].set(5, 4, 7, rand()<<16|rand());//back
		r=sqrt(W2*W2+H2*H2+L2*L2);

		this->density=density;
		set_weight_cuboid(W, L, H);
	}
	void set_weight_cuboid(float W, float L, float H)//dx, dy, dz
	{
		mass=density*W*L*H, inv_mass=1/mass;
		//float m_12=mass/12;
		float m_12=mass/3;
		float W2=W*W, L2=L*L, H2=H*H;
		I.a11=m_12*(L2+H2);//x axis: dy, dz
		I.a22=m_12*(W2+H2);//y axis: dx, dz
		I.a33=m_12*(W2+L2);//z axis: dx, dy
		inv_I=I.inverse_diagonal();
	//	inv_I=I.inverse();
		I_diagonal=true;
	}
};
std::vector<PhysObject> objects0, objects;
int				*txh=0, *txw=0, ntx=0, npols=0, ntpols=0;
int				h=0, w=0, X0, Y0;

Vector3f		cam0(400, 400, 400);
float			dcam0=4;//0.4f
//float const	camx0=400, camy0=400, camz0=400, dcam0=4;
//float const	camx0=4, camy0=4, camz0=4, dcam0=0.04;
float const		ax0=225*_2pi/360, ay0=324.7356103172454f*_2pi/360, tanfov0=1, fov0=.125f*_2pi, da0=2*_2pi/360;
//float const	ax0=3.1532853373782164f, ay0=6.2736568316636330f, tanfov0=0.017029283871059675f, fov0=.125f*_2pi, da0=2*_2pi/360;

Vector3f		cam=cam0;
float			dcam=dcam0;
//float			cam.x=camx0, cam.y=camy0, cam.z=camz0, dcam=dcam0;
float			ax=ax0, ay=ay0, tanfov=tanfov0, fov=fov0, da=da0, da_tfov=tanfov, cax=cos(ax), sax=sin(ax), cay=cos(ay), say=sin(ay);

char			mode=6,//1
				drag=0, d_bypass=0, timer=0, kb[256]={0}, kp=0, line[128], linelen;
tagPOINT		centerP, mouseP0;
tagRECT			R;
HDC__			*ghDC;
HWND__			*ghWnd;
_LARGE_INTEGER	li;
long long		freq, nticks;
int				mod(int x, int y){return x%y-(y&-(x<0));}
inline int		my_abs(int x)
{
	int negative=-!!(0x80000000&x);
	return (x^negative)-negative;
}
inline float	my_abs(float x)
{
	(int&)x&=0x7FFFFFFF;
	return x;
}
//template<typename T>inline T clamp_positive(T x){return 0.5*(x+abs(x));}
//template<int>int clamp_positive(int x){return (x+my_abs(x))>>1;}
inline int		clamp_positive(int x){return (x+my_abs(x))>>1;}
inline float	clamp_positive(float x){return (x+my_abs(x))*0.5f;}
inline float	clamp_negative(float x){return (x-my_abs(x))*0.5f;}
int				minimum(int a, int b){return (a+b-my_abs(a-b))>>1;}
int				minimum(int a, int b, int c)
{
	int a2=a<<1, temp=b+c-my_abs(b-c);
	return (a2+temp-my_abs(a2-temp))>>2;
}
int				maximum(int a, int b){return (a+b+my_abs(a-b))>>1;}
float			clamp(float x, float min, float max)
{
	float t=0.5f*(x+min+abs(x-min));
	return 0.5f*(t+max-abs(t-max));
//	return x<min?min:x>max?max:x;
}
float			interpolate_line_y(float x1, float y1, float x2, float y2, float x){return (y2-y1)/(x2-x1)*(x-x1)+y1;}
float			interpolate_line_x(float x1, float y1, float x2, float y2, float y){return (x2-x1)/(y2-y1)*(y-y1)+x1;}
bool			check_sse4_1()//https://stackoverflow.com/questions/6121792/how-to-check-if-a-cpu-supports-the-sse3-instruction-set
{
	int info[4];
	__cpuid(info, 0);
	int nIds=info[0];
	__cpuid(info, 0x80000000);
	unsigned nExIds=info[0];
	if(nIds>=1)
	{
		__cpuid(info, 1);
		return (info[2]&(1<<19))!=0;
	}
	return false;
}
	const bool	SSE4_1=check_sse4_1();
//	const bool	SSE4_1=false;

//float inv_sqrt(float x)//http://stackoverflow.com/questions/11513344/how-to-implement-the-fast-inverse-square-root-in-java
//{
//	float t0;
//	(long long&)t0=0x5FE6EC85E7DE30DA-((long long&)x>>1);
//	return t0*(1.5-.5*x*t0*t0);
//}

void			world_to_camera_absolute(Vector3f const &p, Vector3f &p_cam)
{
	float cpt=p.x*cax+p.y*sax;
	p_cam.set(p.x*sax-p.y*cax, cpt*say-p.z*cay, cpt*cay+p.z*say);
}
void			world_to_camera(Vector3f const &p, Vector3f &p_cam)
{
	Vector3f d=p-cam;
//	float dx=p.x-cam.x, dy=p.y-cam.y, dz=p.z-cam.z;
	float cpt=d.x*cax+d.y*sax;
	p_cam.set(d.x*sax-d.y*cax, cpt*say-d.z*cay, cpt*cay+d.z*say);
}
void			camera_to_screen(Vector3f const &cp, Vector2f &s)
{
	float temp=cp.z*tanfov/X0;
	s.set(X0+cp.x/temp, Y0+cp.y/temp);
}
void			world_to_screen(Vector3f const &p_world, Vector2f &p_screen, Vector3f &p_cam)
{
	Vector3f d=p_world-cam;
//	float dx=p_world.x-cam.x, dy=p_world.y-cam.y, dz=p_world.z-cam.z;
	float cpt=d.x*cax+d.y*sax;
	p_cam.set(d.x*sax-d.y*cax, cpt*say-d.z*cay, cpt*cay+d.z*say);
	cpt=p_cam.z*tanfov/X0;
	p_screen.set(X0+p_cam.x/cpt, Y0+p_cam.y/cpt);
}
void			world_to_screen(Vector3f const &p, float *v)//*v: Xs, Ys, Xcp, Ycp, Zcp
{
	Vector3f d=p-cam;
//	float dx=p.x-cam.x, dy=p.y-cam.y, dz=p.z-cam.z;
	float cpt=d.x*cax+d.y*sax;
	v[2]=d.x*sax-d.y*cax, v[3]=cpt*say-d.z*cay, v[4]=cpt*cay+d.z*say;
	cpt=v[4]*tanfov/X0;
	v[0]=X0+v[2]/cpt, v[1]=Y0+v[3]/cpt;
}
void			world_to_camera(Vector3f const &p, float *v)//*v: Xcp, Ycp, Zcp
{
	Vector3f d=p-cam;
//	float dx=p.x-cam.x, dy=p.y-cam.y, dz=p.z-cam.z;
	float cpt=d.x*cax+d.y*sax;
	v[0]=d.x*sax-d.y*cax, v[1]=cpt*say-d.z*cay, v[2]=cpt*cay+d.z*say;
}
const int		g_bufsize=1024;
int				g_buflen;
char			g_buf[g_bufsize];
void			GUIPrint(HDC__ *ghDC, int x, int y, const char *a, ...)
{
	g_buflen=vsprintf_s(g_buf, g_bufsize, a, (char*)(&a+1));
//	g_buflen=vsprintf_s(g_buf, g_bufsize, a, (char*)(&reinterpret_cast<const char &>(a))+((sizeof(a)+3)&~3));
	if(g_buflen>0)
		TextOutA(ghDC, x, y, g_buf, g_buflen);
}
void			GUIPrint(HDC__ *ghDC, int x, int y, int value)
{
	g_buflen=sprintf_s(g_buf, 1024, "%d", value);
	if(g_buflen>0)
		TextOutA(ghDC, x, y, g_buf, g_buflen);
}
void			GUIVPrint(HDC__ *ghDC, int x, int y, const char *a, char *va_list)
{
	g_buflen=vsprintf_s(g_buf, g_bufsize, a, va_list);
	if(g_buflen>0)
		TextOutA(ghDC, x, y, g_buf, g_buflen);
}

struct			Mode
{
	virtual void initiate()=0;
	virtual void resize()=0;
	virtual void changeFov()=0;
	virtual void clear_screen()=0;
	virtual void print(int x, int y, const char *a, ...)=0;
	virtual void print(Vector3f const &p, const char *a, ...)=0;
	virtual void draw_point(Vector3f const &p, int color)=0;
	virtual void draw_line(float x1, float y1, float x2, float y2, int lineColor)=0;
	virtual void draw_line(Vector3f const &p1, Vector3f const &p2, int lineColor)=0;
	virtual void draw_triangle(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, int color, bool transparent)=0;
	virtual void draw_ground()=0;
	virtual void render()=0;
	virtual void show()=0;
	virtual void pushTexture(int *texture_)=0;
	virtual int* popTexture()=0;
	virtual void enqueueTextures()=0;
	virtual void clearTextures()=0;
	virtual void finish()=0;
} *m;
struct			LinearSW:private Mode
{
	int **texture;
	int lstart, lend, *libuffer, *lfbuffer;
	float *wbuffer;
//	float *rgb_r, *rgb_g, *rgb_b;
//	int *rgb_r, *rgb_g, *rgb_b;
//	int *rgb_n;
	int *rgb, rgbn;
	HDC__ *ghMemDC;
	HBITMAP__ *hBitmap;
	void initiate();
	void resize();
	void changeFov();
	void clear_screen();
	void print(int x, int y, const char *a, ...);
	void print(Vector3f const &p, const char *a, ...);

	void draw_point(Vector3f const &p, int color);
	void draw_line(float x1, float y1, float x2, float y2, int lineColor);
	void line_A_coeff_x(float Xcp1, float Ycp1, float Zcp1, float Xcp2, float Ycp2, float Zcp2, float &a, float &b);
	void line_A_coeff_y(float Xcp1, float Ycp1, float Zcp1, float Xcp2, float Ycp2, float Zcp2, float &a, float &b);
	void point(float x, float y, float z, int color);
	void _3dLineOnScreen_draw	(Vector2f &s1, Vector3f const &cp1, Vector2f &s2, Vector3f const &cp2, int lineColor);
	void _3dLineOnScreen		(Vector2f &s1, Vector3f const &cp1, Vector2f &s2, Vector3f const &cp2, int lineColor);
	void _2dExtrapolateInfLine	(Vector2f &s1, Vector3f const &cp1, Vector2f &s2, Vector3f const &cp2, int lineColor);
	void _2dExtrapolateInfRay	(Vector2f &s1, Vector3f const &cp1, Vector2f &s2, Vector3f const &cp2, int lineColor);
	void infiniteLine			(Vector3f const &p1, Vector3f const &p2, int lineColor);
	void infiniteRay			(Vector3f const &p1, Vector3f const &p2, int lineColor);
	void _2dExtrapolateLine		(Vector2f &s1, Vector3f const &cp1, Vector2f &s2, Vector3f const &cp2, int lineColor);
	void draw_line				(Vector3f const &p1, Vector3f const &p2, int lineColor);

	void draw_triangle(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, int color, bool transparent);
	void draw_ground();
	void render_opaque		(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, bool textured, int tc, Vector2f const &tx1, Matrix22 const &txm);
	void render_transparent	(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, bool textured, int tc, Vector2f const &tx1, Matrix22 const &txm);
	//void render_opaque(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, int color);
	//void render_transparent(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, int tx_id, Vector2f const &tx1, Matrix22 const &txm);
	void render();
	void show();

	void draft_1behind		(Vector2f &s1, Vector2f &s2, Vector2f &s3);
	void draft_2behind		(Vector2f &s1, Vector2f &s2, Vector2f &s3);
	void draft_start		(Vector2f &s1, Vector2f &s2);
	void draft				(Vector2f &s1, Vector2f &s2);
	void draft_crit_start	(Vector2f &s1, Vector2f &s2);
	void draft_crit			(Vector2f &s1, Vector2f &s2);

	void pushTexture(int *texture_);
	void enqueueTextures();
	int* popTexture();
	void clearTextures();
	void finish();
} lsw;
void				LinearSW::initiate()
{
	texture=0;
	rgbn=h*w;
	wbuffer=(float*)malloc(rgbn*sizeof(float)), libuffer=(int*)malloc(h*sizeof(int)), lfbuffer=(int*)malloc(h*sizeof(int));
	ghMemDC=CreateCompatibleDC(ghDC);
//	rgb_r=(float*)malloc(w*h*sizeof(float)), rgb_g=(float*)malloc(w*h*sizeof(float)), rgb_b=(float*)malloc(w*h*sizeof(float)), rgb_n=(int*)malloc(w*h*sizeof(int));
//	rgb_r=(int*)malloc(w*h*sizeof(int)), rgb_g=(int*)malloc(w*h*sizeof(int)), rgb_b=(int*)malloc(w*h*sizeof(int)), rgb_n=(int*)malloc(w*h*sizeof(int));
	tagBITMAPINFO bmpInfo={{sizeof(tagBITMAPINFOHEADER), w, -h, 1, 32, BI_RGB, 0, 0, 0, 0, 0}};
	hBitmap=(HBITMAP__*)SelectObject(ghMemDC, CreateDIBSection(0, &bmpInfo, DIB_RGB_COLORS, (void**)&rgb, 0, 0));
	SetBkMode(ghMemDC, TRANSPARENT);
}
void				LinearSW::resize()
{
	rgbn=h*w;
	wbuffer=(float*)realloc(wbuffer, rgbn*sizeof(float)), libuffer=(int*)realloc(libuffer, h*sizeof(int)), lfbuffer=(int*)realloc(lfbuffer, h*sizeof(int));
	DeleteObject((HBITMAP__*)SelectObject(ghMemDC, hBitmap));
//	rgb_r=(float*)realloc(rgb_r, w*h*sizeof(float)), rgb_g=(float*)realloc(rgb_g, w*h*sizeof(float)), rgb_b=(float*)realloc(rgb_b, w*h*sizeof(float)), rgb_n=(int*)realloc(rgb_n, w*h*sizeof(int));
//	rgb_r=(int*)realloc(rgb_r, w*h*sizeof(int)), rgb_g=(int*)realloc(rgb_g, w*h*sizeof(int)), rgb_b=(int*)realloc(rgb_b, w*h*sizeof(int)), rgb_n=(int*)realloc(rgb_n, w*h*sizeof(int));
	tagBITMAPINFO bmpInfo={{sizeof(tagBITMAPINFOHEADER), w, -h, 1, 32, BI_RGB, 0, 0, 0, 0, 0}};
	hBitmap=(HBITMAP__*)SelectObject(ghMemDC, CreateDIBSection(0, &bmpInfo, DIB_RGB_COLORS, (void**)&rgb, 0, 0));
	SetBkMode(ghMemDC, TRANSPARENT);
}
void				LinearSW::changeFov(){}
void				LinearSW::clear_screen(){memset(rgb, 0xFF, rgbn*sizeof(int)), memset(wbuffer, 0, rgbn*sizeof(float));}
void				LinearSW::print(int x, int y, const char *a, ...){GUIVPrint(ghMemDC, x, y, a, (char*)(&a+1));}
void				LinearSW::print(Vector3f const &p, const char *a, ...)
{
	Vector2f s;
	Vector3f c;
	world_to_screen(p, s, c);
	if(c.z>0)
		GUIVPrint(ghMemDC, (int)s.x, (int)s.y, a, (char*)(&a+1));
}
void				LinearSW::draw_point(Vector3f const &p, int color){point(p.x, p.y, p.z, color);}
void				LinearSW::draw_line(float x1, float y1, float x2, float y2, int lineColor)
{
	Vector3f cp1(x1, y1, 1), cp2(x2, y2, 1);
	_3dLineOnScreen_draw(Vector2f(x1, y1), cp1, Vector2f(x2, y2), cp2, lineColor);
}
void				LinearSW::line_A_coeff_x(float Xcp1, float Ycp1, float Zcp1, float Xcp2, float Ycp2, float Zcp2, float &a, float &b)
{
	float t=(Xcp2-Xcp1)*Zcp1-Xcp1*(Zcp2-Zcp1);
	a=(Zcp1-Zcp2)*tanfov/(X0*t), b=((Zcp2-Zcp1)*tanfov+Xcp2-Xcp1)/t;
}
void				LinearSW::line_A_coeff_y(float Xcp1, float Ycp1, float Zcp1, float Xcp2, float Ycp2, float Zcp2, float &a, float &b)
{
	float t=(Ycp2-Ycp1)*Zcp1-Ycp1*(Zcp2-Zcp1);
	a=(Zcp1-Zcp2)*tanfov/(X0*t), b=((Zcp2-Zcp1)*Y0*tanfov/X0+Ycp2-Ycp1)/t;
}
void				LinearSW::point(float x, float y, float z, int color)
{
	int bx1=0, bx2=w, by1=0, by2=h, bw=w;
	float dx=x-cam.x, dy=y-cam.y, dz=z-cam.z, cpt=dx*cax+dy*sax, Xcp1=dx*sax-dy*cax, Ycp1=cpt*say-dz*cay;
	float Zcp1=cpt*cay+dz*say;
	if(Zcp1>0)
	{
		float Acp1=1/Zcp1;
		cpt=X0/(Zcp1*tanfov);
		float Xs1=X0+Xcp1*cpt, Ys1=Y0+Ycp1*cpt;
		if((Xs1>=bx1)&(Xs1<bx2)&(Ys1>=by1)&(Ys1<by2))
		{
			int pos=int(Xs1)+bw*int(Ys1);
			if(Acp1>wbuffer[pos])
				rgb[pos]=color, wbuffer[pos]=Acp1;
		}
	}
}
void				LinearSW::_3dLineOnScreen_draw(Vector2f &s1, Vector3f const &cp1, Vector2f &s2, Vector3f const &cp2, int lineColor)
{
	float &x1=s1.x, &y1=s1.y, &x2=s2.x, &y2=s2.y;
	int bx1=0, bx2=w, by1=0, by2=h;
	float xa, ya, xb, yb, dx=x2-x1, dy=y2-y1;
	if(abs(dx)>abs(dy))//horizontal
	{
		float dy_dx=dy/dx;
		if(x1<x2)
		{
			if(x1<bx1)
				xa=(float)bx1, ya=y1+dy_dx*((float)bx1-x1);
			else
				xa=x1, ya=y1;
			if(x2>bx2-1)
				xb=(float)bx2-1, yb=y1+dy_dx*((float)bx2-1-x1);
			else
				xb=x2, yb=y2;
		}
		else
		{
			if(x2<bx1)
				xa=(float)bx1, ya=y1+dy_dx*((float)bx1-x1);
			else
				xa=x2, ya=y2;
			if(x1>bx2-1)
				xb=(float)bx2-1, yb=y1+dy_dx*((float)bx2-1-x1);
			else
				xb=x1, yb=y1;
		}
		float a, b;
		line_A_coeff_x(cp1.x, cp1.y, cp1.z, cp2.x, cp2.y, cp2.z, a, b);
	/*	{
			float t=(Xcp2-Xcp1)*Zcp1-Xcp1*(Zcp2-Zcp1);
			a=(Zcp1-Zcp2)*tanfov/(X0*t), b=((Zcp2-Zcp1)*tanfov+Xcp2-Xcp1)/t;
		}*/
		if(SSE4_1)
		{
			int xEnd=int(xb-xa)-(int(xb-xa)&7);
			for(int x=0;x<xEnd;x+=4)
			{
				__m128 xf=_mm_set_ps(float(x+3), float(x+2), float(x+1), (float)x);
				xf=_mm_add_ps(xf, _mm_set1_ps((float)xa));
				__m128 yf=_mm_sub_ps(xf, _mm_set1_ps((float)x1));
				yf=_mm_mul_ps(yf, _mm_set1_ps((float)dy_dx));
				yf=_mm_add_ps(yf, _mm_set1_ps((float)y1));
				__m128i m_y=_mm_cvtps_epi32(yf);
				__m128i c1=_mm_cmpgt_epi32(m_y, _mm_set1_epi32(by1));//
				__m128i c2=_mm_cmplt_epi32(m_y, _mm_set1_epi32(by2));
				c1=_mm_and_si128(c1, c2);

				c2=_mm_srli_si128(c1, 8);
				c2=_mm_or_si128(c2, c1);
				__m128i c3=_mm_srli_si128(c2, 4);
				c2=_mm_or_si128(c2, c3);
				if(c2.m128i_i32[0])
				{
					__m128i m_x=_mm_cvtps_epi32(xf);
					__m128i pos=_mm_mullo_epi32(m_y, _mm_set1_epi32(w));//SSE4.1
					pos=_mm_add_epi32(pos, m_x);
					__m128 A=_mm_mul_ps(xf, _mm_set1_ps((float)a));
					A=_mm_add_ps(A, _mm_set1_ps((float)b));
				//	__m128 wbk=_mm_set_ps(wbuffer+pos.m128i_i32[0], wbuffer+pos.m128i_i32[0]
					if(c1.m128i_i32[0]&&A.m128_f32[0]>wbuffer[pos.m128i_i32[0]])
						rgb[pos.m128i_i32[0]]=lineColor, wbuffer[pos.m128i_i32[0]]=A.m128_f32[0];
					if(c1.m128i_i32[1]&&A.m128_f32[1]>wbuffer[pos.m128i_i32[1]])
						rgb[pos.m128i_i32[1]]=lineColor, wbuffer[pos.m128i_i32[1]]=A.m128_f32[1];
					if(c1.m128i_i32[2]&&A.m128_f32[2]>wbuffer[pos.m128i_i32[2]])
						rgb[pos.m128i_i32[2]]=lineColor, wbuffer[pos.m128i_i32[2]]=A.m128_f32[2];
					if(c1.m128i_i32[3]&&A.m128_f32[3]>wbuffer[pos.m128i_i32[3]])
						rgb[pos.m128i_i32[3]]=lineColor, wbuffer[pos.m128i_i32[3]]=A.m128_f32[3];
				}
			}
			for(int x=xEnd<0?0:xEnd, xEnd2=int(xb-xa);x<=xEnd2;++x)//horizontal
			{
				int xx=x+(int)xa;
				int y=int(std::floor(y1+dy_dx*(xx-x1)));//-0.5 truncated as 0
				if(y>=by1&&y<by2)
				{
					int pos=w*y+xx;
					float A=a*xx+b;
					if(A>wbuffer[pos])
						rgb[pos]=lineColor, wbuffer[pos]=A;
					else if(A==wbuffer[pos])
					{
						auto p=(unsigned char*)&rgb[pos], c=(unsigned char*)&lineColor;//little endian
						p[0]=(p[0]+c[0])>>1;//b
						p[1]=(p[1]+c[1])>>1;//g
						p[2]=(p[2]+c[2])>>1;//r
					}
				}
			}
		}
		else
		{
			for(int x=int(xa), xEnd=int(xb);x<=xEnd;++x)//horizontal
			{
				//int y;
				//{
				//	float Y=y1+dy_dx*(x-x1);
				//	y=int(Y)-(Y<0);
				//}
				int y=int(std::floor(y1+dy_dx*(x-x1)));//-0.5 truncated as 0
				if(y>=by1&&y<by2)
				{
					int pos=w*y+x;
					float A=a*x+b;
					if(A>wbuffer[pos])
						rgb[pos]=lineColor, wbuffer[pos]=A;
					else if(A==wbuffer[pos])
					{
						auto p=(unsigned char*)&rgb[pos], c=(unsigned char*)&lineColor;//little endian
						p[0]=(p[0]+c[0])>>1;//b
						p[1]=(p[1]+c[1])>>1;//g
						p[2]=(p[2]+c[2])>>1;//r
					}
					//{
					//	((unsigned char*)&rgb[pos])[0]=((unsigned char*)&rgb[pos])[0]+((unsigned char*)&lineColor)[0]>>1;
					//	((unsigned char*)&rgb[pos])[1]=((unsigned char*)&rgb[pos])[1]+((unsigned char*)&lineColor)[1]>>1;
					//	((unsigned char*)&rgb[pos])[2]=((unsigned char*)&rgb[pos])[2]+((unsigned char*)&lineColor)[2]>>1;
					//}
					//	rgb[pos]=0xFFC0CB;
				}
			}
		}
	}
	else//vertical
	{
		float dx_dy=dx/dy;
		if(y1<y2)
		{
			if(y1<by1)
				xa=x1+dx_dy*((float)by1-y1), ya=(float)by1;
			else
				xa=x1, ya=y1;
			if(y2>by2-1)
				xb=x1+dx_dy*((float)by2-1-y1), yb=(float)by2-1;
			else
				xb=x2, yb=y2;
		}
		else
		{
			if(y2<by1)
				xa=x1+dx_dy*((float)by1-y1), ya=(float)by1;
			else
				xa=x2, ya=y2;
			if(y1>by2-1)
				xb=x1+dx_dy*((float)by2-1-y1), yb=(float)by2-1;
			else
				xb=x1, yb=y1;
		}
		float a, b;
		line_A_coeff_y(cp1.x, cp1.y, cp1.z, cp2.x, cp2.y, cp2.z, a, b);
	/*	{
			float t=(Ycp2-Ycp1)*Zcp1-Ycp1*(Zcp2-Zcp1);
			a=(Zcp1-Zcp2)*tanfov/(X0*t), b=((Zcp2-Zcp1)*Y0*tanfov/X0+Ycp2-Ycp1)/t;
		}*/
		if(SSE4_1)
		{
			int yEnd=int(yb-ya)-(int(yb-ya)&7);
			for(int y=0;y<yEnd;y+=4)
			{
				__m128 yf=_mm_set_ps(float(y+3), float(y+2), float(y+1), (float)y);
				yf=_mm_add_ps(yf, _mm_set1_ps((float)ya));
				__m128 xf=_mm_sub_ps(yf, _mm_set1_ps((float)y1));
				xf=_mm_mul_ps(xf, _mm_set1_ps((float)dx_dy));
				xf=_mm_add_ps(xf, _mm_set1_ps((float)x1));
				__m128i m_x=_mm_cvtps_epi32(xf);
				__m128i c1=_mm_cmpgt_epi32(m_x, _mm_set1_epi32(bx1));//
				__m128i c2=_mm_cmplt_epi32(m_x, _mm_set1_epi32(bx2));
				c1=_mm_and_si128(c1, c2);

				c2=_mm_srli_si128(c1, 8);
				c2=_mm_or_si128(c2, c1);
				__m128i c3=_mm_srli_si128(c2, 4);
				c2=_mm_or_si128(c2, c3);
				if(c2.m128i_i32[0])
				{
					__m128i m_y=_mm_cvtps_epi32(yf);
					__m128i pos=_mm_mullo_epi32(m_y, _mm_set1_epi32(w));
					pos=_mm_add_epi32(pos, m_x);
					__m128 A=_mm_mul_ps(yf, _mm_set1_ps((float)a));
					A=_mm_add_ps(A, _mm_set1_ps((float)b));
					if(c1.m128i_i32[0]&&A.m128_f32[0]>wbuffer[pos.m128i_i32[0]])
						rgb[pos.m128i_i32[0]]=lineColor, wbuffer[pos.m128i_i32[0]]=A.m128_f32[0];
					if(c1.m128i_i32[1]&&A.m128_f32[1]>wbuffer[pos.m128i_i32[1]])
						rgb[pos.m128i_i32[1]]=lineColor, wbuffer[pos.m128i_i32[1]]=A.m128_f32[1];
					if(c1.m128i_i32[2]&&A.m128_f32[2]>wbuffer[pos.m128i_i32[2]])
						rgb[pos.m128i_i32[2]]=lineColor, wbuffer[pos.m128i_i32[2]]=A.m128_f32[2];
					if(c1.m128i_i32[3]&&A.m128_f32[3]>wbuffer[pos.m128i_i32[3]])
						rgb[pos.m128i_i32[3]]=lineColor, wbuffer[pos.m128i_i32[3]]=A.m128_f32[3];
				}
			}
			for(int y=yEnd<0?0:yEnd, yEnd2=int(yb-ya);y<=yEnd2;++y)//vertical
			{
				int yy=y+(int)ya;
				int x=int(std::floor(x1+dx_dy*(yy-y1)));//-0.5 truncated as 0
				if(x>=bx1&&x<bx2)
				{
					int pos=w*yy+x;
					float A=a*yy+b;
					if(A>wbuffer[pos])
						rgb[pos]=lineColor, wbuffer[pos]=A;
					else if(A==wbuffer[pos])
					{
						auto p=(unsigned char*)&rgb[pos], c=(unsigned char*)&lineColor;
						p[0]=(p[0]+c[0])>>1;//b
						p[1]=(p[1]+c[1])>>1;//g
						p[2]=(p[2]+c[2])>>1;//r
					}
				}
			}
		}
		else
		{
			for(int y=int(ya), yEnd=int(yb);y<yEnd;++y)//vertical
			{
				int x=int(std::floor(x1+dx_dy*(y-y1)));//-0.5 truncated as 0
				if(x>=bx1&&x<bx2)
				{
					int pos=w*y+x;
					float A=a*y+b;
					if(A>wbuffer[pos])
						rgb[pos]=lineColor, wbuffer[pos]=A;
					else if(A==wbuffer[pos])
					{
						auto p=(unsigned char*)&rgb[pos], c=(unsigned char*)&lineColor;
						p[0]=(p[0]+c[0])>>1;//b
						p[1]=(p[1]+c[1])>>1;//g
						p[2]=(p[2]+c[2])>>1;//r
					}
				}
			}
		}
	}
}
void				LinearSW::_3dLineOnScreen(Vector2f &s1, Vector3f const &cp1, Vector2f &s2, Vector3f const &cp2, int lineColor)
{
	float &x1=s1.x, &y1=s1.y, &x2=s2.x, &y2=s2.y;
	int bx1=0, bx2=w, by1=0, by2=h;
	if(abs(x1)>2e9||abs(y1)>2e9||abs(x2)>2e9||abs(y2)>2e9)
	{
		bool valid=true;
		//up
		if(y1<by1)
		{
			if(y2<by1)
				valid=false;
			else
				x1=interpolate_line_x(x1, y1, x2, y2, (float)by1), y1=(float)by1;
		}
		else if(y2<by1)
			x2=interpolate_line_x(x1, y1, x2, y2, (float)by1), y2=(float)by1;

		if(valid)
		{
			//right
			if(x1>bx2)
			{
				if(x2>bx2)
					valid=false;
				else
					y1=interpolate_line_y(x1, y1, x2, y2, (float)bx2), x1=(float)bx2;
			}
			else if(x2>bx2)
				y2=interpolate_line_y(x1, y1, x2, y2, (float)bx2), x2=(float)bx2;

			if(valid)
			{
				//bottom
				if(y1>by2)
				{
					if(y2>by2)
						valid=false;
					else
						x1=interpolate_line_x(x1, y1, x2, y2, (float)by2), y1=(float)by2;
				}
				else if(y2>by2)
					x2=interpolate_line_x(x1, y1, x2, y2, (float)by2), y2=(float)by2;

				if(valid)
				{
					//left
					if(x1<bx1)
					{
						if(x2<bx1)
							valid=false;
						else
							y1=interpolate_line_y(x1, y1, x2, y2, (float)bx1), x1=(float)bx1;
					}
					else if(x2<bx1)
						y2=interpolate_line_y(x1, y1, x2, y2, (float)bx1), x2=(float)bx1;

					if(valid)
						_3dLineOnScreen_draw(s1, cp1, s2, cp2, lineColor);
				}
			}
		}
	}
	else
		_3dLineOnScreen_draw(s1, cp1, s2, cp2, lineColor);
}
void				LinearSW::_2dExtrapolateInfLine(Vector2f &s1, Vector3f const &cp1, Vector2f &s2, Vector3f const &cp2, int lineColor)
{
	float &x1=s1.x, &y1=s1.y, &x2=s2.x, &y2=s2.y;
	int bx1=0, bx2=w, by1=0, by2=h, bw=w;
	if(x1==x2&&y1==y2)
	{
		if(cp1.z>0)
		{
			float Acp1=1/cp1.z;
			if(x1>=bx1&&x1<bx2&&y1>=by1&&y1<by2)
			{
				int pos=int(x1)+bw*int(y1);
				if(Acp1>wbuffer[pos])
					rgb[pos]=lineColor, wbuffer[pos]=Acp1;
			}
		}
	}
	else
	{
		Vector2f d12(x2-x1, y2-y1);
		float r=sqrt2*(abs(x2)+abs(y2)+2*(X0+Y0))/(abs(d12.x)+abs(d12.y));
		_3dLineOnScreen(s1-r*d12, cp1, s1+r*d12, cp2, lineColor);//[<- 1 2 ->]
	}
}
void				LinearSW::_2dExtrapolateInfRay(Vector2f &s1, Vector3f const &cp1, Vector2f &s2, Vector3f const &cp2, int lineColor)
{
	float &x1=s1.x, &y1=s1.y, &x2=s2.x, &y2=s2.y;
	int bx1=0, bx2=w, by1=0, by2=h, bw=w;
	if(x1==x2&&y1==y2)
	{
		if(cp1.z>0)
		{
			float Acp1=1/cp1.z;
			if(x1>=bx1&&x1<bx2&&y1>=by1&&y1<by2)
			{
				int pos=int(x1)+bw*int(y1);
				if(Acp1>wbuffer[pos])
					rgb[pos]=lineColor, wbuffer[pos]=Acp1;
			}
		}
	}
	else
	{
		Vector2f d12(x2-x1, y2-y1);
		float r=sqrt2*(abs(x2)+abs(y2)+2*(X0+Y0))/(abs(d12.x)+abs(d12.y));
		_3dLineOnScreen(s1, cp1, s2+r*d12, cp2, lineColor);//[1 2 ->]
	}
}
void				LinearSW::infiniteLine(Vector3f const &p1, Vector3f const &p2, int lineColor)//[inf <- 1 2 -> inf]
{
	Vector3f cp1, cpInf;
	world_to_camera(p1, cp1), world_to_camera_absolute(p2-p1, cpInf);//p1, 1->2 from cam
//	float Xcp1, Ycp1, Zcp1;
//	point_world_camera(x1-cam.x, y1-cam.y, z1-cam.z, Xcp1, Ycp1, Zcp1);//P1
//	float XcpInf, YcpInf, ZcpInf;
//	point_world_camera(x2-x1, y2-y1, z2-z1, XcpInf, YcpInf, ZcpInf);//1-2 from cam
	if(cpInf.z>0)//1-2 points forward
	{
		Vector2f sInf;
		camera_to_screen(cpInf, sInf);//Inf
		//float XsInf, YsInf;
		//point_camera_screen(XcpInf, YcpInf, ZcpInf, XsInf, YsInf);//Inf

		Vector3f cp2;
		world_to_camera(p2, cp2);
		//float Xcp2, Ycp2, Zcp2;
		//point_world_camera(x2-cam.x, y2-cam.y, z2-cam.z, Xcp2, Ycp2, Zcp2);//P2
		if(cp1.z>0)//P1 in front of the camera
		{
			Vector2f s1;
			camera_to_screen(cp1, s1);
			//float Xs1, Ys1;
			//point_camera_screen(Xcp1, Ycp1, Zcp1, Xs1, Ys1);
			_2dExtrapolateInfRay(sInf, cp1, s1, cp2, lineColor);
		}
		else
		{
			if(cp1.z==0)//P1 at the camera plane
				cp1+=cp1-cp2;
			//	Xcp1+=Xcp1-Xcp2, Ycp1+=Ycp1-Ycp2, Zcp1+=Zcp1-Zcp2;
			Vector2f s1;
			camera_to_screen(cp1, s1);
			//float Xs1, Ys1;
			//point_camera_screen(Xcp1, Ycp1, Zcp1, Xs1, Ys1);
			_2dExtrapolateLine(s1, cp1, sInf, cp2, lineColor);//P1 behind the camera
		}
	}
	else if(cpInf.z==0)//1-2 is parallel to cam plane
	{
		if(cp1.z>0)//in front of cam
		{
			Vector3f cp2;
			Vector2f s1, s2;
			world_to_screen(p2, s2, cp2);
			camera_to_screen(cp1, s1);
			//float Xcp2, Ycp2, Zcp2, Xs2, Ys2;
			//point_world_screen(x2-cam.x, y2-cam.y, z2-cam.z, Xcp2, Ycp2, Zcp2, Xs2, Ys2);//P2
			//float Xs1, Ys1;
			//point_camera_screen(Xcp1, Ycp1, Zcp1, Xs1, Ys1);
			_2dExtrapolateInfLine(s2, cp1, s1, cp2, lineColor);
		}
	}
	else//1-2 points back
	{
		Vector2f sInf;
		Vector3f cp2;
		camera_to_screen(cpInf, sInf);
		world_to_camera(p2, cp2);
		//float XsInf, YsInf;
		//point_camera_screen(XcpInf, YcpInf, ZcpInf, XsInf, YsInf);//Inf

		//float Xcp2, Ycp2, Zcp2;
		//point_world_camera(x2-cam.x, y2-cam.y, z2-cam.z, Xcp2, Ycp2, Zcp2);//P2
		if(cp1.z>0)//P1 in front of cam
		{
			Vector2f s1;
			camera_to_screen(cp1, s1);
			//float Xs1, Ys1;
			//point_camera_screen(Xcp1, Ycp1, Zcp1, Xs1, Ys1);
			_2dExtrapolateInfRay(sInf, cp1, s1, cp2, lineColor);
		}
		else
		{
			if(cp1.z==0)//P1 at cam plane
				cp1-=cp1-cp2;
			Vector2f s1;
			camera_to_screen(cp1, s1);
			//if(Zcp1==0)
			//	Xcp1-=Xcp1-Xcp2, Ycp1-=Ycp1-Ycp2, Zcp1-=Zcp1-Zcp2;
			//float Xs1, Ys1;
			//point_camera_screen(Xcp1, Ycp1, Zcp1, Xs1, Ys1);
			_2dExtrapolateLine(s1, cp1, sInf, cp2, lineColor);//P1 behind the camera
		}
	}
}
void				LinearSW::infiniteRay(Vector3f const &p1, Vector3f const &p2, int lineColor)//[1 2 -> inf]
{
	Vector3f cp1, cpInf;
	world_to_camera(p1, cp1), world_to_camera_absolute(p2-p1, cpInf);//P1: starting point, Inf: direction from camera
	//float Xcp1, Ycp1, Zcp1;
	//point_world_camera(x1-cam.x, y1-cam.y, z1-cam.z, Xcp1, Ycp1, Zcp1);
	//float XcpInf, YcpInf, ZcpInf;
	//point_world_camera(x2-x1, y2-y1, z2-z1, XcpInf, YcpInf, ZcpInf);//
	if(cpInf.z>0)//ray pointing forward
	{
		Vector2f sInf;
		camera_to_screen(cpInf, sInf);//Inf
		Vector3f cp2;
		world_to_camera(p2, cp2);//P2
		//float XsInf, YsInf;
		//point_camera_screen(XcpInf, YcpInf, ZcpInf, XsInf, YsInf);
		//float Xcp2, Ycp2, Zcp2;
		//point_world_camera(x2-cam.x, y2-cam.y, z2-cam.z, Xcp2, Ycp2, Zcp2);
		if(cp1.z>0)//starts in front of the camera
		{
			Vector2f s1;
			camera_to_screen(cp1, s1);
			//float Xs1, Ys1;
			//point_camera_screen(Xcp1, Ycp1, Zcp1, Xs1, Ys1);
			_3dLineOnScreen(s1, cp1, sInf, cp2, lineColor);
		}
		else if(cp1.z==0)//starts at the camera plane
		{
			cp1+=cp1-cp2;
			Vector2f s1;
			camera_to_screen(cp1, s1);//1 behind
			//Xcp1+=Xcp1-Xcp2, Ycp1+=Ycp1-Ycp2, Zcp1+=Zcp1-Zcp2;
			//float Xs1, Ys1;
			//point_camera_screen(Xcp1, Ycp1, Zcp1, Xs1, Ys1);
			_2dExtrapolateLine(s1, cp1, sInf, cp2, lineColor);
		}
		else//starts behind the camera
		{
			Vector2f s1;
			camera_to_screen(cp1, s1);
			//float Xs1, Ys1;
			//point_camera_screen(Xcp1, Ycp1, Zcp1, Xs1, Ys1);
			_2dExtrapolateLine(s1, cp1, sInf, cp2, lineColor);
		}
	}
	else if(cpInf.z==0)//ray pointing parallel to camera plane
	{
		if(cp1.z>0)//starts in front of the camera
		{
			Vector3f cp2;
			Vector2f s2, s1;
			world_to_screen(p1, s2, cp2);//P2
			//float Xcp2, Ycp2, Zcp2, Xs2, Ys2;
			//point_world_screen(x2-cam.x, y2-cam.y, z2-cam.z, Xcp2, Ycp2, Zcp2, Xs2, Ys2);
			camera_to_screen(cp1, s1);
			//float Xs1, Ys1;
			//point_camera_screen(Xcp1, Ycp1, Zcp1, Xs1, Ys1);
			_2dExtrapolateInfRay(s2, cp1, s1, cp2, lineColor);
		}
	}
	else//ray pointing back
	{
		if(cp1.z>0)//starts in front of the camera
		{
			Vector2f sInf, s1;
			Vector3f cp2;
			camera_to_screen(cpInf, sInf);
			world_to_camera(p2, cp2);
			camera_to_screen(cp1, s1);
			//float XsInf, YsInf;
			//point_camera_screen(XcpInf, YcpInf, ZcpInf, XsInf, YsInf);//Inf
			//float Xcp2, Ycp2, Zcp2;
			//point_world_camera(x2-cam.x, y2-cam.y, z2-cam.z, Xcp2, Ycp2, Zcp2);//P2
			//float Xs1, Ys1;
			//point_camera_screen(Xcp1, Ycp1, Zcp1, Xs1, Ys1);
			_2dExtrapolateLine(sInf, cp1, s1, cp2, lineColor);
		}
	}
}
void				LinearSW::_2dExtrapolateLine(Vector2f &s1, Vector3f const &cp1, Vector2f &s2, Vector3f const &cp2, int lineColor)//P1 behind
{
	float &x1=s1.x, &y1=s1.y, &x2=s2.x, &y2=s2.y;
	int bx1=0, bx2=w, by1=0, by2=h, bw=w;
	if(x1==x2&&y1==y2)
	{
		if(cp1.z>0)
		{
			float Acp1=1/cp1.z;
			if(x1>=bx1&&x1<bx2&&y1>=by1&&y1<by2)
			{
				int pos=int(x1)+bw*int(y1);
				if(Acp1>wbuffer[pos])
					rgb[pos]=lineColor, wbuffer[pos]=Acp1;
			}
		}
	}
	else
	{
		float dx=x2-x1, dy=y2-y1, r=sqrt2*(abs(x2)+abs(y2)+2*(X0+Y0))/(abs(dx)+abs(dy));
		_3dLineOnScreen(s2, cp1, Vector2f(x2+r*dx, y2+r*dy), cp2, lineColor);
	}
}
void				LinearSW::draw_line(Vector3f const &p1, Vector3f const &p2, int lineColor)
{
	int drawn=0;
	const float &x1=p1.x, &y1=p1.y, &z1=p1.z, &x2=p2.x, &y2=p2.y, &z2=p2.z;
	switch(is_nan_or_inf(x1)+is_nan_or_inf(y1)+is_nan_or_inf(z1)+is_nan_or_inf(x2)+is_nan_or_inf(y2)+is_nan_or_inf(z2))
	{
	case 0:
		{
			Vector2f s1, s2;
			Vector3f cp1, cp2;
			world_to_camera(p1, cp1), world_to_camera(p2, cp2);
			//{
			//	Vector3f d=p1-cam;
			//	cp1.set(d.x*sax-d.y*cax, cpt*say-d.z*cay, cpt*cay+d.z*say)
			////	float dx=x1-cam.x, dy=y1-cam.y, dz=z1-cam.z, cpt=dx*cax+dy*sax;
			////	Xcp1=dx*sax-dy*cax, Ycp1=cpt*say-dz*cay, Zcp1=cpt*cay+dz*say;
			//}
		//	float dx=x2-cam.x, dy=y2-cam.y, dz=z2-cam.z, cpt=dx*cax+dy*sax, Xcp2=dx*sax-dy*cax, Ycp2=cpt*say-dz*cay, Zcp2=cpt*cay+dz*say;
			if(cp2.z>0)
			{
				if(cp1.z>0)
				{
					camera_to_screen(cp1, s1), camera_to_screen(cp2, s2);
					//float	temp=X0/(cp1.z*tanfov);	s1.x=X0+cp1.x*temp, s1.y=Y0+cp1.y*temp;
					//		temp=X0/(cp2.z*tanfov);	s2.x=X0+cp2.x*temp, s2.y=Y0+cp2.y*temp;
					//cpt=X0/(Zcp1*tanfov), Xs1=X0+Xcp1*cpt, Ys1=Y0+Ycp1*cpt;
					//cpt=X0/(Zcp2*tanfov), Xs2=X0+Xcp2*cpt, Ys2=Y0+Ycp2*cpt;
					_3dLineOnScreen(s1, cp1, s2, cp2, lineColor);
				}
				else if(cp1.z==0)
				{
					cp1+=cp1-cp2;
					camera_to_screen(cp1, s1), camera_to_screen(cp2, s2);
					//Xcp1+=Xcp1-Xcp2, Ycp1+=Ycp1-Ycp2, Zcp1+=Zcp1-Zcp2;
					//cpt=X0/(Zcp1*tanfov), Xs1=X0+Xcp1*cpt, Ys1=Y0+Ycp1*cpt;
					//cpt=X0/(Zcp2*tanfov), Xs2=X0+Xcp2*cpt, Ys2=Y0+Ycp2*cpt;
					_2dExtrapolateLine(s1, cp1, s2, cp2, lineColor);
				}
				else
				{
					camera_to_screen(cp1, s1), camera_to_screen(cp2, s2);
					//cpt=X0/(Zcp1*tanfov), Xs1=X0+Xcp1*cpt, Ys1=Y0+Ycp1*cpt;
					//cpt=X0/(Zcp2*tanfov), Xs2=X0+Xcp2*cpt, Ys2=Y0+Ycp2*cpt;
					_2dExtrapolateLine(s1, cp1, s2, cp2, lineColor);
				}
				drawn=1;
			}
			else if(cp2.z==0)
			{
				if(cp1.z>0)
				{
					cp2+=cp2-cp1;
					camera_to_screen(cp1, s1), camera_to_screen(cp2, s2);
					//Xcp2+=Xcp2-Xcp1, Ycp2+=Ycp2-Ycp1, Zcp2+=Zcp2-Zcp1;
					//cpt=X0/(Zcp1*tanfov), Xs1=X0+Xcp1*cpt, Ys1=Y0+Ycp1*cpt;
					//cpt=X0/(Zcp2*tanfov), Xs2=X0+Xcp2*cpt, Ys2=Y0+Ycp2*cpt;
					_2dExtrapolateLine(s2, cp1, s1, cp2, lineColor), drawn=1;
				}
			}
			else if(cp1.z>0)
			{
				camera_to_screen(cp1, s1), camera_to_screen(cp2, s2);
				//cpt=X0/(Zcp1*tanfov), Xs1=X0+Xcp1*cpt, Ys1=Y0+Ycp1*cpt;
				//cpt=X0/(Zcp2*tanfov), Xs2=X0+Xcp2*cpt, Ys2=Y0+Ycp2*cpt;
				_2dExtrapolateLine(s2, cp1, s1, cp2, lineColor), drawn=1;
			}
		}
		break;
	case 1:
		{
			if(isinf(x1))
				infiniteRay(p2, Vector3f(x1==infinity?x2+1:x2-1, y2, z2), lineColor), drawn=1;
			else if(isinf(y1))
				infiniteRay(p2, Vector3f(x2, y1==infinity?y2+1:y2-1, z2), lineColor), drawn=1;
			else if(isinf(z1))
				infiniteRay(p2, Vector3f(x2, y2, z1==infinity?z2+1:z2-1), lineColor), drawn=1;
			else if(isinf(x2))
				infiniteRay(p1, Vector3f(x2==infinity?x1+1:x1-1, y1, z1), lineColor), drawn=1;
			else if(isinf(y2))
				infiniteRay(p1, Vector3f(x1, y2==infinity?y1+1:y1-1, z1), lineColor), drawn=1;
			else if(isinf(z2))
				infiniteRay(p1, Vector3f(x1, y1, z2==infinity?z1+1:z1-1), lineColor), drawn=1;
		}
		break;
	case 2:
		{
			if(isinf(x1)&&isinf(x2))
			{
				if(y1==x2&&z1==z2)
					infiniteLine(Vector3f(cam.x, x2, z2), Vector3f(cam.x+1, x2, z2), lineColor), drawn=1;
			}
			else if(isinf(y1)&&isinf(x2))
			{
				if(z1==z2&&x1==x2)
					infiniteLine(Vector3f(x2, cam.y, z2), Vector3f(x2, cam.y+1, z2), lineColor), drawn=1;
			}
			else if(isinf(z1)&&isinf(z2))
			{
				if(x1==x2&&y1==x2)
					infiniteLine(Vector3f(x2, x2, cam.z), Vector3f(x2, x2, cam.z+1), lineColor), drawn=1;
			}
		}
		break;
	}
	if(!drawn)
	{
		if(!(is_nan_or_inf(x1)+is_nan_or_inf(y1)+is_nan_or_inf(z1)))
			point(x1, y1, z1, lineColor);
		else if(!(is_nan_or_inf(x2)+is_nan_or_inf(x2)+is_nan_or_inf(z2)))
			point(x1, y1, z1, lineColor);
	}
}
void				LinearSW::draw_triangle(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, int color, bool transparent)
{
	if(transparent)
		render_transparent(p1, p2, p3, false, color, Vector2f(), Matrix22());
	else
		render_opaque(p1, p2, p3, false, color, Vector2f(), Matrix22());
}
void				LinearSW::draw_ground()
{
	if(!ntx)
	{
		txh=(int*)realloc(txh, sizeof(int)), txh[ntx]=100;
		txw=(int*)realloc(txw, sizeof(int)), txw[ntx]=100;
		int *texture_=(int*)malloc(txh[ntx]*txw[ntx]*sizeof(int));
		for(int ky=0, kyEnd=txh[ntx];ky<kyEnd;++ky)
			for(int kx=0, kxEnd=txw[ntx];kx<kxEnd;++kx)
				texture_[kxEnd*ky+kx]=(unsigned char)(127+127*sin(kx*_2pi/kxEnd))<<16|(unsigned char)(127+127*sin(ky*_2pi/kyEnd))<<8|(unsigned char)(127+127*sin(ky*_2pi/kyEnd*kx*_2pi/kxEnd));
		//for(int k=0;k<txh[ntx]*txw[ntx];++k)
		//	texture_[k]=(rand()<<15|rand())&0xFFFFFF;
		pushTexture(texture_);
	//	enqueueTextures();
	}
	Vector3f h_arrow(cam.x+cax, cam.y+sax, cam.z), cp1;//vector to horizon
	Vector2f s1;
	world_to_screen(h_arrow, s1, cp1);

	//Vector2f s1, s2, s3;
	Vector3f c1, c2, c3;
	float X0_tfov=X0/tanfov;
	world_to_camera(Vector3f(0, 0, 0), c1), world_to_camera(Vector3f(1000, 0, 0), c2), world_to_camera(Vector3f(0, 1000, 0), c3);
	memset(libuffer, 0, h*sizeof(int));
	for(int k=0;k<h;++k)lfbuffer[k]=w;
//	world_to_screen(p1, s1, c1), world_to_screen(p2, s2, c2), world_to_screen(p3, s3, c3);
	//		if(c1.z<0){		 if(c2.z<0){		 if(c3.z<0)	return;
	//										else			draft_2behind(s3, s1, s2);}
	//					else{					 if(c3.z<0)	draft_2behind(s2, s3, s1);
	//										else			draft_1behind(s2, s3, s1);}}
	//else{					 if(c2.z<0){		 if(c3.z<0)	draft_2behind(s1, s2, s3);
	//										else			draft_1behind(s3, s1, s2);}
	//					else{					 if(c3.z<0)	draft_1behind(s1, s2, s3);
	//										else			draft_start(s1, s2), draft(s2, s3), draft(s3, s1);}}
	Vector3f upy=c2-c1,		//u12	=<12>
		upx=c3-c1,			//ux3	=<13>
		n=upy.cross(upx);	//abc	=<n>	=<12>x<13>
	float t=n.dot(c1);
	if(!t)//camera is on the plane
		return;
	float B7=n.x, B8=n.y, B9=n.z*X0_tfov-n.x*X0-n.y*Y0;	float cpt=t*X0_tfov;
	float A1=B7/cpt, A2=B8/cpt, A3=B9/cpt;

	float	ux3_1=upy.magnitude();	upy/=ux3_1;			//u12	=<u12>	=<12>/|12|
			ux3_1=upx.dot(upy),		upx-=ux3_1*upy;		//ux3	=<x3>	=<13>-<13>.<u12><u12>
			ux3_1=upx.magnitude(),	upx/=ux3_1;			//ux3	=<ux3>	=<x3>/|x3|
			ux3_1=upx.dot(c1);							//ux3_1			=<ux3>.<1>
	float	u12_1=upy.dot(c1);							//u12_1			=<u12>.<1>

	float tx_multiplier=0.01f;
	Vector2f C1=tx_multiplier*Vector2f(upx.x, upy.x), C2=tx_multiplier*Vector2f(upx.y, upy.y), C3=tx_multiplier*Vector2f(upx.z, upy.z);
	float	D1=C1.x*t, D2=C2.x*t, D3=(C3.x*X0_tfov-C1.x*X0-C2.x*Y0)*t;
	float	D4=C1.y*t, D5=C2.y*t, D6=(C3.y*X0_tfov-C1.y*X0-C2.y*Y0)*t;
	Vector2f E1=-tx_multiplier*Vector2f(ux3_1, u12_1);
	float	B1=D1+E1.x*B7, B2=D2+E1.x*B8, B3=D3+E1.x*B9;
	float	B4=D4+E1.y*B7, B5=D5+E1.y*B8, B6=D6+E1.y*B9;
	float *wb_y=wbuffer;
	int y_start=int(s1.y>=h?h:s1.y<0?0:s1.y);
	int lib_y, lfb_y, *rgb_y=rgb+y_start*w, *txk=texture[0], txhk=txh[0], txwk=txw[0];
	if(SSE4_1)
	{
		__m128 adm_increment=_mm_set1_ps(A1*4);
		__m128 m_a_increment=_mm_set1_ps(B1*4);
		__m128 m_b_increment=_mm_set1_ps(B4*4);
		__m128 m_c_increment=_mm_set1_ps(B7*4);

		__m128 s_txwk=_mm_set1_ps((float)txwk);
		__m128 s_txhk=_mm_set1_ps((float)txhk);
		__m128i i_txwk=_mm_set1_epi32(txwk), i_txwk_1=_mm_set1_epi32(txwk-1);
		__m128i i_txhk=_mm_set1_epi32(txhk), i_txhk_1=_mm_set1_epi32(txhk-1);
		for(int ys=y_start;ys<h;++ys)
		{
			lib_y=libuffer[ys]<0?0:libuffer[ys], lfb_y=lfbuffer[ys]>w?w:lfbuffer[ys];
			float admittance=A1*lib_y+A2*ys+A3, a=B1*lib_y+B2*ys+B3, b=B4*lib_y+B5*ys+B6, c=B7*lib_y+B8*ys+B9;

			int xs=lib_y;
			__m128 adm=_mm_set_ps(admittance+3*A1, admittance+A1+A1, admittance+A1, admittance);
			__m128 m_a=_mm_set_ps(a+B1*3, a+B1+B1, a+B1, a);
			__m128 m_b=_mm_set_ps(b+B4*3, b+B4+B4, b+B4, b);
			__m128 m_c=_mm_set_ps(c+B7*3, c+B7+B7, c+B7, c);
			for(int xsEnd=lfb_y-4;xs<xsEnd;xs+=4)//for each 4 pixels
			{
				__m128 wb=_mm_loadu_ps(wb_y+xs);
				__m128i c_result=_mm_castps_si128(_mm_cmplt_ps(wb, adm));
				__m128i c_result2=_mm_shuffle_epi32(c_result, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
				__m128i c_or=_mm_or_si128(c_result, c_result2);//r3|r1, r2|r0, r1|r3, r0|r2
				c_result2=_mm_shuffle_epi32(c_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
				c_or=_mm_or_si128(c_or, c_result2);//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
				if(c_or.m128i_i32[0])
				{
					__m128 wb2=_mm_castsi128_ps(_mm_and_si128(c_result, _mm_castps_si128(adm)));
					c_result2=_mm_xor_si128(c_result, _mm_set1_epi32(-1));
					wb=_mm_castsi128_ps(_mm_and_si128(_mm_castps_si128(wb), c_result2));
					wb2=_mm_castsi128_ps(_mm_or_si128(_mm_castps_si128(wb2), _mm_castps_si128(wb)));
					_mm_storeu_ps(wb_y+xs, wb2);
					__m128 a_c=_mm_div_ps(m_a, m_c);
					__m128 b_c=_mm_div_ps(m_b, m_c);

					//mod operation		x%q = x-q*floor(x/q)
					__m128 ac_w=_mm_div_ps(a_c, s_txwk);
					__m128 bc_h=_mm_div_ps(b_c, s_txhk);

					__m128 f_acw=_mm_floor_ps(ac_w);//floor(x/txw)	SSE4.1
					__m128 f_bch=_mm_floor_ps(bc_h);

					f_acw=_mm_mul_ps(f_acw, s_txwk);//floor(x/txw)*txw
					f_bch=_mm_mul_ps(f_bch, s_txhk);

					a_c=_mm_sub_ps(a_c, f_acw);//x-floor(x/txw)*txw
					b_c=_mm_sub_ps(b_c, f_bch);
						
					a_c=_mm_floor_ps(a_c);//99.9 -> 99 not 100	SSE4.1
					b_c=_mm_floor_ps(b_c);

					__m128i xtx=_mm_cvtps_epi32(a_c);
					__m128i ytx=_mm_cvtps_epi32(b_c);

					__m128i cmp_x=_mm_cmplt_epi32(xtx, _mm_setzero_si128());//index guard
					__m128i cmp_y=_mm_cmplt_epi32(ytx, _mm_setzero_si128());
					cmp_x=_mm_and_si128(cmp_x, i_txwk);
					cmp_y=_mm_and_si128(cmp_y, i_txhk);
					xtx=_mm_add_epi32(xtx, cmp_x);
					ytx=_mm_add_epi32(ytx, cmp_y);
					cmp_x=_mm_cmpgt_epi32(xtx, i_txwk_1);
					cmp_y=_mm_cmpgt_epi32(ytx, i_txhk_1);
					cmp_x=_mm_and_si128(cmp_x, i_txwk);
					cmp_y=_mm_and_si128(cmp_y, i_txhk);
					xtx=_mm_sub_epi32(xtx, cmp_x);
					ytx=_mm_sub_epi32(ytx, cmp_y);

					__m128i tx_idx=_mm_mullo_epi32(i_txwk, ytx);//SSE4.1
					tx_idx=_mm_add_epi32(tx_idx, xtx);

					int txs=txwk*txhk;
					if(tx_idx.m128i_i32[3]<0||tx_idx.m128i_i32[3]>=txs||tx_idx.m128i_i32[2]<0||tx_idx.m128i_i32[2]>=txs||tx_idx.m128i_i32[1]<0||tx_idx.m128i_i32[1]>=txs||tx_idx.m128i_i32[0]<0||tx_idx.m128i_i32[0]>=txs)
						int LOL_1=0;
					__m128i m_tx=_mm_set_epi32(txk[tx_idx.m128i_i32[3]], txk[tx_idx.m128i_i32[2]], txk[tx_idx.m128i_i32[1]], txk[tx_idx.m128i_i32[0]]);
					__m128i m_rgb=_mm_loadu_si128((__m128i*)(rgb_y+xs));
					m_tx=_mm_and_si128(m_tx, c_result);
					m_rgb=_mm_and_si128(m_rgb, c_result2);
					m_tx=_mm_or_si128(m_tx, m_rgb);
					_mm_storeu_si128((__m128i*)(rgb_y+xs), m_tx);
				}
				adm=_mm_add_ps(adm, adm_increment);
				m_a=_mm_add_ps(m_a, m_a_increment);
				m_b=_mm_add_ps(m_b, m_b_increment);
				m_c=_mm_add_ps(m_c, m_c_increment);
			}
			admittance=adm.m128_f32[0];
			a=m_a.m128_f32[0];
			b=m_b.m128_f32[0];
			c=m_c.m128_f32[0];
			for(;xs<lfb_y;++xs)
			{
				if(admittance>wb_y[xs])
				{
					wb_y[xs]=admittance;
					int Xtx=int(a/c)%txwk;Xtx+=txwk&-(Xtx<0);
					int Ytx=int(b/c)%txhk;Ytx+=txhk&-(Ytx<0);
					rgb_y[xs]=txk[Xtx+Ytx*txwk];
				}
				admittance+=A1, a+=B1, b+=B4, c+=B7;
			}
			wb_y=wb_y+w, rgb_y=rgb_y+w;
		}
	}
	else//no SSE4.1
	{
		const __m128 half=_mm_set1_ps(0.5f);
		const __m128i mask_lo=_mm_set_epi32(0, -1, 0, -1), mask_hi=_mm_set_epi32(-1, 0, -1, 0);
		__m128 s_txwk=_mm_set1_ps((float)txwk);
		__m128 s_txhk=_mm_set1_ps((float)txhk);
		__m128i i_txwk=_mm_set1_epi32(txwk), i_txwk_1=_mm_set1_epi32(txwk-1);
		__m128i i_txhk=_mm_set1_epi32(txhk), i_txhk_1=_mm_set1_epi32(txhk-1);

		__m128 adm_increment=_mm_set1_ps(A1*4);
		__m128 m_a_increment=_mm_set1_ps(B1*4);
		__m128 m_b_increment=_mm_set1_ps(B4*4);
		__m128 m_c_increment=_mm_set1_ps(B7*4);
		for(int ys=y_start;ys<h;++ys)
		{
			lib_y=libuffer[ys]<0?0:libuffer[ys], lfb_y=lfbuffer[ys]>w?w:lfbuffer[ys];
			float admittance=A1*lib_y+A2*ys+A3, a=B1*lib_y+B2*ys+B3, b=B4*lib_y+B5*ys+B6, c=B7*lib_y+B8*ys+B9;

			int xs=lib_y;
			__m128 adm=_mm_set_ps(admittance+3*A1, admittance+A1+A1, admittance+A1, admittance);
			__m128 m_a=_mm_set_ps(a+B1*3, a+B1+B1, a+B1, a);
			__m128 m_b=_mm_set_ps(b+B4*3, b+B4+B4, b+B4, b);
			__m128 m_c=_mm_set_ps(c+B7*3, c+B7+B7, c+B7, c);
			for(int xsEnd=lfb_y-4;xs<xsEnd;xs+=4)//for each 4 pixels
			{
				__m128 wb=_mm_loadu_ps(wb_y+xs);
				__m128i c_result=_mm_castps_si128(_mm_cmplt_ps(wb, adm));
				__m128i c_result2=_mm_shuffle_epi32(c_result, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
				__m128i c_or=_mm_or_si128(c_result, c_result2);//r3|r1, r2|r0, r1|r3, r0|r2
				c_result2=_mm_shuffle_epi32(c_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
				c_or=_mm_or_si128(c_or, c_result2);//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
				if(c_or.m128i_i32[0])
				{
					__m128 wb2=_mm_castsi128_ps(_mm_and_si128(c_result, _mm_castps_si128(adm)));
					c_result2=_mm_xor_si128(c_result, _mm_set1_epi32(-1));
					wb=_mm_castsi128_ps(_mm_and_si128(_mm_castps_si128(wb), c_result2));
					wb2=_mm_castsi128_ps(_mm_or_si128(_mm_castps_si128(wb2), _mm_castps_si128(wb)));
					_mm_storeu_ps(wb_y+xs, wb2);
					__m128 a_c=_mm_div_ps(m_a, m_c);
					__m128 b_c=_mm_div_ps(m_b, m_c);

					//mod operation		x%q = x-q*floor(x/q)
					__m128 ac_w=_mm_div_ps(a_c, s_txwk);
					__m128 bc_h=_mm_div_ps(b_c, s_txhk);

					//__m128 f_acw=_mm_floor_ps(ac_w);//floor(x/txw)	SSE4.1
					//__m128 f_bch=_mm_floor_ps(bc_h);
					ac_w=_mm_sub_ps(ac_w, half);//floor(x/q)
					bc_h=_mm_sub_ps(bc_h, half);
					__m128i i_acw=_mm_cvtps_epi32(ac_w);
					__m128i i_bch=_mm_cvtps_epi32(bc_h);
					__m128 f_acw=_mm_cvtepi32_ps(i_acw);
					__m128 f_bch=_mm_cvtepi32_ps(i_bch);

					f_acw=_mm_mul_ps(f_acw, s_txwk);//floor(x/txw)*txw
					f_bch=_mm_mul_ps(f_bch, s_txhk);

					a_c=_mm_sub_ps(a_c, f_acw);//x-floor(x/txw)*txw
					b_c=_mm_sub_ps(b_c, f_bch);
						
					//a_c=_mm_floor_ps(a_c);//99.9 -> 99 not 100	SSE4.1
					//b_c=_mm_floor_ps(b_c);
					a_c=_mm_sub_ps(a_c, half);
					b_c=_mm_sub_ps(b_c, half);

					__m128i xtx=_mm_cvtps_epi32(a_c);
					__m128i ytx=_mm_cvtps_epi32(b_c);

					__m128i cmp_x=_mm_cmplt_epi32(xtx, _mm_setzero_si128());//index guard
					__m128i cmp_y=_mm_cmplt_epi32(ytx, _mm_setzero_si128());
					cmp_x=_mm_and_si128(cmp_x, i_txwk);
					cmp_y=_mm_and_si128(cmp_y, i_txhk);
					xtx=_mm_add_epi32(xtx, cmp_x);
					ytx=_mm_add_epi32(ytx, cmp_y);
					cmp_x=_mm_cmpgt_epi32(xtx, i_txwk_1);
					cmp_y=_mm_cmpgt_epi32(ytx, i_txhk_1);
					cmp_x=_mm_and_si128(cmp_x, i_txwk);
					cmp_y=_mm_and_si128(cmp_y, i_txhk);
					xtx=_mm_sub_epi32(xtx, cmp_x);
					ytx=_mm_sub_epi32(ytx, cmp_y);

					//__m128i tx_idx=_mm_mullo_epi32(i_txwk, ytx);//SSE4.1
					__m128i tx_idx_lo=_mm_mul_epu32(ytx, i_txwk);
					ytx=_mm_shuffle_epi32(ytx, _MM_SHUFFLE(2, 3, 0, 1));//y2 y3 y0 y1
					__m128i tx_idx_hi=_mm_mul_epu32(ytx, i_txwk);
					tx_idx_hi=_mm_shuffle_epi32(tx_idx_hi, _MM_SHUFFLE(2, 3, 0, 1));
					tx_idx_lo=_mm_and_si128(tx_idx_lo, mask_lo);
					tx_idx_hi=_mm_and_si128(tx_idx_hi, mask_hi);
					__m128i tx_idx=_mm_or_si128(tx_idx_lo, tx_idx_hi);

					tx_idx=_mm_add_epi32(tx_idx, xtx);

					int txs=txwk*txhk;
					if(tx_idx.m128i_i32[3]<0||tx_idx.m128i_i32[3]>=txs||tx_idx.m128i_i32[2]<0||tx_idx.m128i_i32[2]>=txs||tx_idx.m128i_i32[1]<0||tx_idx.m128i_i32[1]>=txs||tx_idx.m128i_i32[0]<0||tx_idx.m128i_i32[0]>=txs)
						int LOL_1=0;
					__m128i m_tx=_mm_set_epi32(txk[tx_idx.m128i_i32[3]], txk[tx_idx.m128i_i32[2]], txk[tx_idx.m128i_i32[1]], txk[tx_idx.m128i_i32[0]]);
					__m128i m_rgb=_mm_loadu_si128((__m128i*)(rgb_y+xs));
					m_tx=_mm_and_si128(m_tx, c_result);
					m_rgb=_mm_and_si128(m_rgb, c_result2);
					m_tx=_mm_or_si128(m_tx, m_rgb);
					_mm_storeu_si128((__m128i*)(rgb_y+xs), m_tx);
				}
				adm=_mm_add_ps(adm, adm_increment);
				m_a=_mm_add_ps(m_a, m_a_increment);
				m_b=_mm_add_ps(m_b, m_b_increment);
				m_c=_mm_add_ps(m_c, m_c_increment);
			}
			admittance=adm.m128_f32[0];
			a=m_a.m128_f32[0];
			b=m_b.m128_f32[0];
			c=m_c.m128_f32[0];
			for(;xs<lfb_y;++xs)
			{
				if(admittance>wb_y[xs])
				{
					wb_y[xs]=admittance;
					int Xtx=int(a/c)%txwk;Xtx+=txwk&-(Xtx<0);
					int Ytx=int(b/c)%txhk;Ytx+=txhk&-(Ytx<0);
					rgb_y[xs]=txk[Xtx+Ytx*txwk];
				}
				admittance+=A1, a+=B1, b+=B4, c+=B7;
			}
			wb_y=wb_y+w, rgb_y=rgb_y+w;
		}
		//for(int ys=0;ys<h;++ys)
		//{
		//	lib_y=libuffer[ys]<0?0:libuffer[ys], lfb_y=lfbuffer[ys]>w?w:lfbuffer[ys];
		//	float admittance=A1*lib_y+A2*ys+A3, a=B1*lib_y+B2*ys+B3, b=B4*lib_y+B5*ys+B6, c=B7*lib_y+B8*ys+B9;
		//	for(int xs=lib_y;xs<lfb_y;++xs)
		//	{
		//		if(std::abs(admittance-wb_y[xs])<1e-10*std::abs(admittance+wb_y[xs]))//z fighting
		//		{
		//			int Xtx=int(a/c)%txwk;Xtx+=txwk&-(Xtx<0);
		//			int Ytx=int(b/c)%txhk;Ytx+=txhk&-(Ytx<0);
		//			auto ps=(unsigned char*)&rgb_y[xs], pt=(unsigned char*)&txk[Xtx+Ytx*txwk];
		//			ps[0]=(ps[0]+pt[0])>>1;//b	warning C4554
		//			ps[1]=(ps[1]+pt[1])>>1;//g
		//			ps[2]=(ps[2]+pt[2])>>1;//r
		//		}
		//		else if(admittance>wb_y[xs])
		//		{
		//			wb_y[xs]=admittance;
		//			int Xtx=int(a/c)%txwk;Xtx+=txwk&-(Xtx<0);
		//			int Ytx=int(b/c)%txhk;Ytx+=txhk&-(Ytx<0);
		//			rgb_y[xs]=txk[Xtx+Ytx*txwk];
		//		}
		//		admittance+=A1, a+=B1, b+=B4, c+=B7;
		//	}
		//	wb_y=wb_y+w, rgb_y=rgb_y+w;
		//}
	}
}
void				LinearSW::render_opaque		(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, bool textured, int tc, Vector2f const &tx1, Matrix22 const &txm)
{
	Vector2f s1, s2, s3;
	Vector3f c1, c2, c3;
	float X0_tfov=X0/tanfov;
	world_to_screen(p1, s1, c1), world_to_screen(p2, s2, c2), world_to_screen(p3, s3, c3);
			if(c1.z<0){		 if(c2.z<0){		 if(c3.z<0)	return;
											else			draft_2behind(s3, s1, s2);}
						else{					 if(c3.z<0)	draft_2behind(s2, s3, s1);
											else			draft_1behind(s2, s3, s1);}}
	else{					 if(c2.z<0){		 if(c3.z<0)	draft_2behind(s1, s2, s3);
											else			draft_1behind(s3, s1, s2);}
						else{					 if(c3.z<0)	draft_1behind(s1, s2, s3);
											else			draft_start(s1, s2), draft(s2, s3), draft(s3, s1);}}
	Vector3f upy=c2-c1,		//u12	=<12>
		upx=c3-c1,			//ux3	=<13>
		a=upy.cross(upx);	//abc	=<n>	=<12>x<13>
	float t=a.dot(c1);
	if(!t)
		return;
	float B7=a.x, B8=a.y, B9=a.z*X0_tfov-a.x*X0-a.y*Y0;	float cpt=t*X0_tfov;
	float A1=B7/cpt, A2=B8/cpt, A3=B9/cpt;
	if(textured)//textured triangle
	{
		float	ux3_1=upy.magnitude();	upy/=ux3_1;			//u12	=<u12>	=<12>/|12|
				ux3_1=upx.dot(upy),		upx-=ux3_1*upy;		//ux3	=<x3>	=<13>-<13>.<u12><u12>
				ux3_1=upx.magnitude(),	upx/=ux3_1;			//ux3	=<ux3>	=<x3>/|x3|
				ux3_1=upx.dot(c1);							//ux3_1			=<ux3>.<1>
		float	u12_1=upy.dot(c1);							//u12_1			=<u12>.<1>

		Vector2f C1=txm*Vector2f(upx.x, upy.x), C2=txm*Vector2f(upx.y, upy.y), C3=txm*Vector2f(upx.z, upy.z);
		float	D1=C1.x*t, D2=C2.x*t, D3=(C3.x*X0_tfov-C1.x*X0-C2.x*Y0)*t;
		float	D4=C1.y*t, D5=C2.y*t, D6=(C3.y*X0_tfov-C1.y*X0-C2.y*Y0)*t;
		Vector2f E1=tx1-txm*Vector2f(ux3_1, u12_1);
		float	B1=D1+E1.x*B7, B2=D2+E1.x*B8, B3=D3+E1.x*B9;
		float	B4=D4+E1.y*B7, B5=D5+E1.y*B8, B6=D6+E1.y*B9;
		float *wb_y=wbuffer;
		int lib_y, lfb_y, *rgb_y=rgb, *txk=texture[tc], txhk=txh[tc], txwk=txw[tc];
		if(SSE4_1)
		{
			const __m128 half=_mm_set1_ps(0.5f);
			const __m128i mask_lo=_mm_set_epi32(0, 0, -1, -1), mask_hi=_mm_set_epi32(-1, -1, 0, 0);
			__m128i i_txwk_1=_mm_set1_epi32(txwk-1);
			__m128i i_txhk_1=_mm_set1_epi32(txhk-1);

			__m128 s_txwk=_mm_set1_ps((float)txwk);
			__m128 s_txhk=_mm_set1_ps((float)txhk);
			__m128i i_txwk=_mm_set1_epi32(txwk);
			__m128i i_txhk=_mm_set1_epi32(txhk);
			__m128 adm_increment=_mm_set1_ps(A1*4);
			__m128 m_a_increment=_mm_set1_ps(B1*4);
			__m128 m_b_increment=_mm_set1_ps(B4*4);
			__m128 m_c_increment=_mm_set1_ps(B7*4);
			for(int ys=0;ys<h;++ys)
			{
				lib_y=libuffer[ys]<0?0:libuffer[ys], lfb_y=lfbuffer[ys]>w?w:lfbuffer[ys];
				float admittance=A1*lib_y+A2*ys+A3, a=B1*lib_y+B2*ys+B3, b=B4*lib_y+B5*ys+B6, c=B7*lib_y+B8*ys+B9;

				int xs=lib_y;
				__m128 adm=_mm_set_ps(admittance+3*A1, admittance+A1+A1, admittance+A1, admittance);
				__m128 m_a=_mm_set_ps(a+B1*3, a+B1+B1, a+B1, a);
				__m128 m_b=_mm_set_ps(b+B4*3, b+B4+B4, b+B4, b);
				__m128 m_c=_mm_set_ps(c+B7*3, c+B7+B7, c+B7, c);
				for(int xsEnd=lfb_y-4;xs<xsEnd;xs+=4)//for each 4 pixels
				{
					__m128 wb=_mm_loadu_ps(wb_y+xs);
					__m128i c_result=_mm_castps_si128(_mm_cmplt_ps(wb, adm));
					__m128i c_result2=_mm_shuffle_epi32(c_result, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
					__m128i c_or=_mm_or_si128(c_result, c_result2);//r3|r1, r2|r0, r1|r3, r0|r2
					c_result2=_mm_shuffle_epi32(c_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
					c_or=_mm_or_si128(c_or, c_result2);//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
					if(c_or.m128i_i32[0])
					{
						__m128 wb2=_mm_castsi128_ps(_mm_and_si128(c_result, _mm_castps_si128(adm)));
						c_result2=_mm_xor_si128(c_result, _mm_set1_epi32(-1));
						wb=_mm_castsi128_ps(_mm_and_si128(_mm_castps_si128(wb), c_result2));
						wb2=_mm_castsi128_ps(_mm_or_si128(_mm_castps_si128(wb2), _mm_castps_si128(wb)));
						_mm_storeu_ps(wb_y+xs, wb2);
						__m128 a_c=_mm_div_ps(m_a, m_c);
						__m128 b_c=_mm_div_ps(m_b, m_c);

						//mod operation		x%q = x-q*floor(x/q)
						__m128 ac_w=_mm_div_ps(a_c, s_txwk);
						__m128 bc_h=_mm_div_ps(b_c, s_txhk);

						__m128 f_acw=_mm_floor_ps(ac_w);//floor(x/txw)	SSE4.1
						__m128 f_bch=_mm_floor_ps(bc_h);

						f_acw=_mm_mul_ps(f_acw, s_txwk);//floor(x/txw)*txw
						f_bch=_mm_mul_ps(f_bch, s_txhk);

						a_c=_mm_sub_ps(a_c, f_acw);//x-floor(x/txw)*txw
						b_c=_mm_sub_ps(b_c, f_bch);
						
						a_c=_mm_floor_ps(a_c);//99.9 -> 99 not 100	SSE4.1
						b_c=_mm_floor_ps(b_c);
						__m128i xtx=_mm_cvtps_epi32(a_c);
						__m128i ytx=_mm_cvtps_epi32(b_c);

						__m128i cmp_x=_mm_cmplt_epi32(xtx, _mm_setzero_si128());//index guard
						__m128i cmp_y=_mm_cmplt_epi32(ytx, _mm_setzero_si128());
						cmp_x=_mm_and_si128(cmp_x, i_txwk);
						cmp_y=_mm_and_si128(cmp_y, i_txhk);
						xtx=_mm_add_epi32(xtx, cmp_x);
						ytx=_mm_add_epi32(ytx, cmp_y);
						cmp_x=_mm_cmpgt_epi32(xtx, i_txwk_1);
						cmp_y=_mm_cmpgt_epi32(ytx, i_txhk_1);
						cmp_x=_mm_and_si128(cmp_x, i_txwk);
						cmp_y=_mm_and_si128(cmp_y, i_txhk);
						xtx=_mm_sub_epi32(xtx, cmp_x);
						ytx=_mm_sub_epi32(ytx, cmp_y);

						__m128i tx_idx=_mm_mullo_epi32(i_txwk, ytx);//SSE4.1
						tx_idx=_mm_add_epi32(tx_idx, xtx);

						__m128i m_tx=_mm_set_epi32(txk[tx_idx.m128i_i32[3]], txk[tx_idx.m128i_i32[2]], txk[tx_idx.m128i_i32[1]], txk[tx_idx.m128i_i32[0]]);
						__m128i m_rgb=_mm_loadu_si128((__m128i*)(rgb_y+xs));
						m_tx=_mm_and_si128(m_tx, c_result);
						m_rgb=_mm_and_si128(m_rgb, c_result2);
						m_tx=_mm_or_si128(m_tx, m_rgb);
						_mm_storeu_si128((__m128i*)(rgb_y+xs), m_tx);
					}
					adm=_mm_add_ps(adm, adm_increment);
					m_a=_mm_add_ps(m_a, m_a_increment);
					m_b=_mm_add_ps(m_b, m_b_increment);
					m_c=_mm_add_ps(m_c, m_c_increment);
				}
				admittance=adm.m128_f32[0];
				a=m_a.m128_f32[0];
				b=m_b.m128_f32[0];
				c=m_c.m128_f32[0];
				for(;xs<lfb_y;++xs)
				{
					//if(std::abs(admittance-wb_y[xs])<1e-10*std::abs(admittance+wb_y[xs]))//z fighting
					//{
					//	int Xtx=int(a/c)%txwk;Xtx+=txwk&-(Xtx<0);//if(Xtx<0)Xtx+=txwk;
					//	int Ytx=int(b/c)%txhk;Ytx+=txhk&-(Ytx<0);//if(Ytx<0)Ytx+=txhk;
					//	auto ps=(unsigned char*)&rgb_y[xs], pt=(unsigned char*)&txk[Xtx+Ytx*txwk];
					//	ps[0]=(ps[0]+pt[0])>>1;//b	warning C4554
					//	ps[1]=(ps[1]+pt[1])>>1;//g
					//	ps[2]=(ps[2]+pt[2])>>1;//r
					//}
					//else
					if(admittance>wb_y[xs])
					{
						wb_y[xs]=admittance;
						int Xtx=int(a/c)%txwk;Xtx+=txwk&-(Xtx<0);//if(Xtx<0)Xtx+=txwk;
						int Ytx=int(b/c)%txhk;Ytx+=txhk&-(Ytx<0);//if(Ytx<0)Ytx+=txhk;
						rgb_y[xs]=txk[Xtx+Ytx*txwk];
					}
					admittance+=A1, a+=B1, b+=B4, c+=B7;
				}
				wb_y=wb_y+w, rgb_y=rgb_y+w;
			}
		}
		else//textured, no SSE4.1
		{
			const __m128 half=_mm_set1_ps(0.5f);
			const __m128i mask2_lo=_mm_set_epi32(0, -1, 0, -1), mask2_hi=_mm_set_epi32(-1, 0, -1, 0);
			__m128 s_txwk=_mm_set1_ps((float)txwk);
			__m128 s_txhk=_mm_set1_ps((float)txhk);
			__m128i i_txwk=_mm_set1_epi32(txwk), i_txwk_1=_mm_set1_epi32(txwk-1);
			__m128i i_txhk=_mm_set1_epi32(txhk), i_txhk_1=_mm_set1_epi32(txhk-1);
			__m128 adm_increment=_mm_set1_ps(A1*4);
			__m128 m_a_increment=_mm_set1_ps(B1*4);
			__m128 m_b_increment=_mm_set1_ps(B4*4);
			__m128 m_c_increment=_mm_set1_ps(B7*4);
			for(int ys=0;ys<h;++ys)
			{
				lib_y=libuffer[ys]<0?0:libuffer[ys], lfb_y=lfbuffer[ys]>w?w:lfbuffer[ys];
				float admittance=A1*lib_y+A2*ys+A3, a=B1*lib_y+B2*ys+B3, b=B4*lib_y+B5*ys+B6, c=B7*lib_y+B8*ys+B9;

				int xs=lib_y;
				__m128 adm=_mm_set_ps(admittance+3*A1, admittance+A1+A1, admittance+A1, admittance);
				__m128 m_a=_mm_set_ps(a+B1*3, a+B1+B1, a+B1, a);
				__m128 m_b=_mm_set_ps(b+B4*3, b+B4+B4, b+B4, b);
				__m128 m_c=_mm_set_ps(c+B7*3, c+B7+B7, c+B7, c);
				for(int xsEnd=lfb_y-4;xs<xsEnd;xs+=4)//for each 4 pixels
				{
					__m128 wb=_mm_loadu_ps(wb_y+xs);
					__m128i c_result=_mm_castps_si128(_mm_cmplt_ps(wb, adm));
					__m128i c_result2=_mm_shuffle_epi32(c_result, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
					__m128i c_or=_mm_or_si128(c_result, c_result2);//r3|r1, r2|r0, r1|r3, r0|r2
					c_result2=_mm_shuffle_epi32(c_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
					c_or=_mm_or_si128(c_or, c_result2);//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
					if(c_or.m128i_i32[0])
					{
						__m128 wb2=_mm_castsi128_ps(_mm_and_si128(c_result, _mm_castps_si128(adm)));
						c_result2=_mm_xor_si128(c_result, _mm_set1_epi32(-1));
						wb=_mm_castsi128_ps(_mm_and_si128(_mm_castps_si128(wb), c_result2));
						wb2=_mm_castsi128_ps(_mm_or_si128(_mm_castps_si128(wb2), _mm_castps_si128(wb)));
						_mm_storeu_ps(wb_y+xs, wb2);
						__m128 a_c=_mm_div_ps(m_a, m_c);
						__m128 b_c=_mm_div_ps(m_b, m_c);

						//mod operation		x%q = x-q*floor(x/q)
						__m128 ac_w=_mm_div_ps(a_c, s_txwk);//x/q
						__m128 bc_h=_mm_div_ps(b_c, s_txhk);

						ac_w=_mm_sub_ps(ac_w, half);//floor(x/q)
						bc_h=_mm_sub_ps(bc_h, half);
						__m128i i_acw=_mm_cvtps_epi32(ac_w);
						__m128i i_bch=_mm_cvtps_epi32(bc_h);
						ac_w=_mm_cvtepi32_ps(i_acw);
						bc_h=_mm_cvtepi32_ps(i_bch);

						ac_w=_mm_mul_ps(ac_w, s_txwk);//floor(x/q)*q
						bc_h=_mm_mul_ps(bc_h, s_txhk);

						a_c=_mm_sub_ps(a_c, ac_w);//x-floor(x/q)*q
						b_c=_mm_sub_ps(b_c, bc_h);

						a_c=_mm_sub_ps(a_c, half);
						b_c=_mm_sub_ps(b_c, half);
						__m128i xtx=_mm_cvtps_epi32(a_c);
						__m128i ytx=_mm_cvtps_epi32(b_c);

						__m128i cmp_x=_mm_cmplt_epi32(xtx, _mm_setzero_si128());//index guard
						__m128i cmp_y=_mm_cmplt_epi32(ytx, _mm_setzero_si128());
						cmp_x=_mm_and_si128(cmp_x, i_txwk);
						cmp_y=_mm_and_si128(cmp_y, i_txhk);
						xtx=_mm_add_epi32(xtx, cmp_x);
						ytx=_mm_add_epi32(ytx, cmp_y);
						cmp_x=_mm_cmpgt_epi32(xtx, i_txwk_1);
						cmp_y=_mm_cmpgt_epi32(ytx, i_txhk_1);
						cmp_x=_mm_and_si128(cmp_x, i_txwk);
						cmp_y=_mm_and_si128(cmp_y, i_txhk);
						xtx=_mm_sub_epi32(xtx, cmp_x);
						ytx=_mm_sub_epi32(ytx, cmp_y);
						
						//__m128i tx_idx=_mm_mullo_epi32(i_txwk, ytx);//SSE4.1
						__m128i tx_idx_lo=_mm_mul_epu32(ytx, i_txwk);
						ytx=_mm_shuffle_epi32(ytx, _MM_SHUFFLE(2, 3, 0, 1));//y2 y3 y0 y1
						__m128i tx_idx_hi=_mm_mul_epu32(ytx, i_txwk);
						tx_idx_hi=_mm_shuffle_epi32(tx_idx_hi, _MM_SHUFFLE(2, 3, 0, 1));
						tx_idx_lo=_mm_and_si128(tx_idx_lo, mask2_lo);
						tx_idx_hi=_mm_and_si128(tx_idx_hi, mask2_hi);
						__m128i tx_idx=_mm_or_si128(tx_idx_lo, tx_idx_hi);

						tx_idx=_mm_add_epi32(tx_idx, xtx);

						//int txs=txwk*txhk;
						//if(tx_idx.m128i_i32[3]<0||tx_idx.m128i_i32[3]>=txs||tx_idx.m128i_i32[2]<0||tx_idx.m128i_i32[2]>=txs||tx_idx.m128i_i32[1]<0||tx_idx.m128i_i32[1]>=txs||tx_idx.m128i_i32[0]<0||tx_idx.m128i_i32[0]>=txs)
						//	int LOL_1=0;
						__m128i m_tx=_mm_set_epi32(txk[tx_idx.m128i_i32[3]], txk[tx_idx.m128i_i32[2]], txk[tx_idx.m128i_i32[1]], txk[tx_idx.m128i_i32[0]]);
						__m128i m_rgb=_mm_loadu_si128((__m128i*)(rgb_y+xs));
						m_tx=_mm_and_si128(m_tx, c_result);
						m_rgb=_mm_and_si128(m_rgb, c_result2);
						m_tx=_mm_or_si128(m_tx, m_rgb);
						_mm_storeu_si128((__m128i*)(rgb_y+xs), m_tx);
						//rgb_y[xs]=m_tx.m128i_i32[0], rgb_y[xs+1]=m_tx.m128i_i32[1], rgb_y[xs+2]=m_tx.m128i_i32[2], rgb_y[xs+3]=m_tx.m128i_i32[3];
					}
					adm=_mm_add_ps(adm, adm_increment);
					m_a=_mm_add_ps(m_a, m_a_increment);
					m_b=_mm_add_ps(m_b, m_b_increment);
					m_c=_mm_add_ps(m_c, m_c_increment);
				}
				admittance=adm.m128_f32[0];
				a=m_a.m128_f32[0];
				b=m_b.m128_f32[0];
				c=m_c.m128_f32[0];
				for(;xs<lfb_y;++xs)
				{
					if(admittance>wb_y[xs])
					{
						wb_y[xs]=admittance;
						int Xtx=int(a/c)%txwk;Xtx+=txwk&-(Xtx<0);//if(Xtx<0)Xtx+=txwk;
						int Ytx=int(b/c)%txhk;Ytx+=txhk&-(Ytx<0);//if(Ytx<0)Ytx+=txhk;
						rgb_y[xs]=txk[Xtx+Ytx*txwk];
					}
					admittance+=A1, a+=B1, b+=B4, c+=B7;
				}
				wb_y=wb_y+w, rgb_y=rgb_y+w;
			}
		/*	for(int ys=0;ys<h;++ys)
			{
				lib_y=libuffer[ys]<0?0:libuffer[ys], lfb_y=lfbuffer[ys]>w?w:lfbuffer[ys];
				float admittance=A1*lib_y+A2*ys+A3, a=B1*lib_y+B2*ys+B3, b=B4*lib_y+B5*ys+B6, c=B7*lib_y+B8*ys+B9;
				for(int xs=lib_y;xs<lfb_y;++xs)
				{
					if(std::abs(admittance-wb_y[xs])<1e-10*std::abs(admittance+wb_y[xs]))//z fighting
					{
						int Xtx=int(a/c)%txwk;Xtx+=txwk&-(Xtx<0);//if(Xtx<0)Xtx+=txwk;
						int Ytx=int(b/c)%txhk;Ytx+=txhk&-(Ytx<0);//if(Ytx<0)Ytx+=txhk;
						auto ps=(unsigned char*)&rgb_y[xs], pt=(unsigned char*)&txk[Xtx+Ytx*txwk];
						ps[0]=(ps[0]+pt[0])>>1;//b	warning C4554
						ps[1]=(ps[1]+pt[1])>>1;//g
						ps[2]=(ps[2]+pt[2])>>1;//r
					}
					else if(admittance>wb_y[xs])
					{
						wb_y[xs]=admittance;
						int Xtx=int(a/c)%txwk;Xtx+=txwk&-(Xtx<0);//if(Xtx<0)Xtx+=txwk;
						int Ytx=int(b/c)%txhk;Ytx+=txhk&-(Ytx<0);//if(Ytx<0)Ytx+=txhk;
						rgb_y[xs]=txk[Xtx+Ytx*txwk];
					}
					admittance+=A1, a+=B1, b+=B4, c+=B7;
				}
				wb_y=wb_y+w, rgb_y=rgb_y+w;
			}//*/
		}
	}
	else//solid color triangle
	{
		float *wb_y=wbuffer;
		int lib_y, lfb_y, *rgb_y=rgb, color=tc;
	//	if(SSE4_1)
		{
			__m128i m_color=_mm_set1_epi32(color);
			for(int ys=0;ys<h;++ys)
			{
				lib_y=libuffer[ys]<0?0:libuffer[ys], lfb_y=lfbuffer[ys]>w?w:lfbuffer[ys];
				float admittance=A1*lib_y+A2*ys+A3;

				int xs=lib_y;
				__m128 adm=_mm_set_ps(admittance+3*A1, admittance+A1+A1, admittance+A1, admittance);
				__m128 adm_increment=_mm_set1_ps(A1*4);
				for(int xsEnd=lfb_y-4;xs<xsEnd;xs+=4)
				{
					__m128 wb=_mm_loadu_ps(wb_y+xs);
					__m128i c_result=_mm_castps_si128(_mm_cmplt_ps(wb, adm));
					__m128i c_result2=_mm_shuffle_epi32(c_result, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
					__m128i c_or=_mm_or_si128(c_result, c_result2);//r3|r1, r2|r0, r1|r3, r0|r2
					c_result2=_mm_shuffle_epi32(c_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
					c_or=_mm_or_si128(c_or, c_result2);//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
					if(c_or.m128i_i32[0])
					{
						__m128 wb2=_mm_castsi128_ps(_mm_and_si128(c_result, _mm_castps_si128(adm)));
						c_result2=_mm_xor_si128(c_result, _mm_set1_epi32(-1));
						wb=_mm_castsi128_ps(_mm_and_si128(_mm_castps_si128(wb), c_result2));
						wb2=_mm_castsi128_ps(_mm_or_si128(_mm_castps_si128(wb2), _mm_castps_si128(wb)));
						_mm_storeu_ps(wb_y+xs, wb2);

						__m128i m_rgb=_mm_loadu_si128((__m128i*)(rgb_y+xs));
						__m128i m_color2=_mm_and_si128(m_color, c_result);
						m_rgb=_mm_and_si128(m_rgb, c_result2);
						m_color2=_mm_or_si128(m_color2, m_rgb);
						_mm_storeu_si128((__m128i*)(rgb_y+xs), m_color2);
					}
					adm=_mm_add_ps(adm, adm_increment);
				}
				admittance=adm.m128_f32[0];
				for(;xs<lfb_y;++xs)
				{
					//if(std::abs(admittance-wb_y[xs])<1e-10*std::abs(admittance+wb_y[xs]))//z fighting
					//{
					//	auto p=(unsigned char*)&rgb_y[xs], c=(unsigned char*)&color;
					//	p[0]=(p[0]+c[0])>>1;//b		warning C4554
					//	p[1]=(p[1]+c[1])>>1;//g
					//	p[2]=(p[2]+c[2])>>1;//r
					//}
					//else
					if(admittance>wb_y[xs])
						wb_y[xs]=admittance, rgb_y[xs]=color;
					admittance+=A1;
				}
				wb_y=wb_y+w, rgb_y=rgb_y+w;
			}
		}
	/*	else//solid color triangle - no SSE4.1
		{
			for(int ys=0;ys<h;++ys)
			{
				lib_y=libuffer[ys]<0?0:libuffer[ys], lfb_y=lfbuffer[ys]>w?w:lfbuffer[ys];
				float admittance=A1*lib_y+A2*ys+A3;
				for(int xs=lib_y;xs<lfb_y;++xs)
				{
					if(std::abs(admittance-wb_y[xs])<1e-10*std::abs(admittance+wb_y[xs]))//z fighting
					{
						auto p=(unsigned char*)&rgb_y[xs], c=(unsigned char*)&color;
						p[0]=(p[0]+c[0])>>1;//b		warning C4554
						p[1]=(p[1]+c[1])>>1;//g
						p[2]=(p[2]+c[2])>>1;//r
					}
					else if(admittance>wb_y[xs])
						wb_y[xs]=admittance, rgb_y[xs]=color;
					admittance+=A1;
				}
				wb_y=wb_y+w, rgb_y=rgb_y+w;
			}
		}//*/
	}
//	GUIPrint(ghMemDC, (int)v1[0], (int)v1[1], k);
}
void				LinearSW::render_transparent(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, bool textured, int tc, Vector2f const &tx1, Matrix22 const &txm)
{
	Vector2f s1, s2, s3;
	Vector3f c1, c2, c3;
	float X0_tfov=X0/tanfov;
	world_to_screen(p1, s1, c1), world_to_screen(p2, s2, c2), world_to_screen(p3, s3, c3);
			if(c1.z<0){		 if(c2.z<0){		 if(c3.z<0)	return;
											else			draft_2behind(s3, s1, s2);}
						else{					 if(c3.z<0)	draft_2behind(s2, s3, s1);
											else			draft_1behind(s2, s3, s1);}}
	else{					 if(c2.z<0){		 if(c3.z<0)	draft_2behind(s1, s2, s3);
											else			draft_1behind(s3, s1, s2);}
						else{					 if(c3.z<0)	draft_1behind(s1, s2, s3);
											else			draft_start(s1, s2), draft(s2, s3), draft(s3, s1);}}
	Vector3f upy=c2-c1,		//u12	=<12>
		upx=c3-c1,			//ux3	=<13>
		a=upy.cross(upx);	//abc	=<n>	=<12>x<13>
	float t=a.dot(c1);
	if(!t)return;
	float B7=a.x, B8=a.y, B9=a.z*X0_tfov-a.x*X0-a.y*Y0;	float cpt=t*X0_tfov;
	float A1=B7/cpt, A2=B8/cpt, A3=B9/cpt;
	if(textured)//textured transparent triangle
	{
		float	ux3_1=upy.magnitude();	upy/=ux3_1;			//u12	=<u12>	=<12>/|12|
				ux3_1=upx.dot(upy),		upx-=ux3_1*upy;		//ux3	=<x3>	=<13>-<13>.<u12><u12>
				ux3_1=upx.magnitude(),	upx/=ux3_1;			//ux3	=<ux3>	=<x3>/|x3|
				ux3_1=upx.dot(c1);							//ux3_1			=<ux3>.<1>
		float	u12_1=upy.dot(c1);							//u12_1			=<u12>.<1>
		Vector2f C1=txm*Vector2f(upx.x, upy.x), C2=txm*Vector2f(upx.y, upy.y), C3=txm*Vector2f(upx.z, upy.z);
		float	D1=C1.x*t, D2=C2.x*t, D3=(C3.x*X0_tfov-C1.x*X0-C2.x*Y0)*t;
		float	D4=C1.y*t, D5=C2.y*t, D6=(C3.y*X0_tfov-C1.y*X0-C2.y*Y0)*t;
		Vector2f E1=tx1-txm*Vector2f(ux3_1, u12_1);
		float	B1=D1+E1.x*B7, B2=D2+E1.x*B8, B3=D3+E1.x*B9;
		float	B4=D4+E1.y*B7, B5=D5+E1.y*B8, B6=D6+E1.y*B9;
		float *wb_y=wbuffer;
		int lib_y, lfb_y, *rgb_y=rgb, *txk=texture[tc], txhk=txh[tc], txwk=txw[tc];
		if(SSE4_1)
		{
			__m128i const sh_mul_lo=_mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 15, 13, 11, 9, 7, 5, 3, 1);
			__m128i const sh_mul_hi=_mm_set_epi8(15, 13, 11, 9, 7, 5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0);
			const __m128i sh_add_lo=_mm_set_epi8(1, 1, 1, 1, 1, 1, 1, 1, 14, 12, 10, 8, 6, 4, 2, 0);
			const __m128i sh_add_hi=_mm_set_epi8(14, 12, 10, 8, 6, 4, 2, 0, 1, 1, 1, 1, 1, 1, 1, 1);
			__m128i const mask_lo=_mm_set_epi32(0, 0, -1, -1);
			__m128i const mask_hi=_mm_set_epi32(-1, -1, 0, 0);
			for(int ys=0;ys<h;++ys)
			{
				lib_y=libuffer[ys]<0?0:libuffer[ys], lfb_y=lfbuffer[ys]>w?w:lfbuffer[ys];
				float admittance=A1*lib_y+A2*ys+A3, a=B1*lib_y+B2*ys+B3, b=B4*lib_y+B5*ys+B6, c=B7*lib_y+B8*ys+B9;

				int xs=lib_y;
				__m128 adm=_mm_set_ps(admittance+3*A1, admittance+A1+A1, admittance+A1, admittance);
				__m128 adm_increment=_mm_set1_ps(A1*4);
				__m128 m_a=_mm_set_ps(a+B1*3, a+B1+B1, a+B1, a);
				__m128 m_a_increment=_mm_set1_ps(B1*4);
				__m128 m_b=_mm_set_ps(b+B4*3, b+B4+B4, b+B4, b);
				__m128 m_b_increment=_mm_set1_ps(B4*4);
				__m128 m_c=_mm_set_ps(c+B7*3, c+B7+B7, c+B7, c);
				__m128 m_c_increment=_mm_set1_ps(B7*4);

				__m128 s_txwk=_mm_set1_ps((float)txwk);
				__m128 s_txhk=_mm_set1_ps((float)txhk);
				__m128i i_txwk=_mm_set1_epi32(txwk);
				__m128i i_txhk=_mm_set1_epi32(txhk);
				for(int xsEnd=lfb_y-4;xs<xsEnd;xs+=4)//for each 4 pixels
				{
					__m128 wb=_mm_loadu_ps(wb_y+xs);
					__m128i c_result=_mm_castps_si128(_mm_cmplt_ps(wb, adm));
					__m128i c_result2=_mm_shuffle_epi32(c_result, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
					__m128i c_or=_mm_or_si128(c_result, c_result2);//r3|r1, r2|r0, r1|r3, r0|r2
					c_result2=_mm_shuffle_epi32(c_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
								c_or=_mm_or_si128(c_or, c_result2);//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
					if(c_or.m128i_i32[0])
					{
						c_result2=_mm_xor_si128(c_result, _mm_set1_epi32(-1));
						__m128 a_c=_mm_div_ps(m_a, m_c);
						__m128 b_c=_mm_div_ps(m_b, m_c);

						//mod operation		x%q = x-q*floor(x/q)
						__m128 ac_w=_mm_div_ps(a_c, s_txwk);
						__m128 bc_h=_mm_div_ps(b_c, s_txhk);

						__m128 f_acw=_mm_floor_ps(ac_w);//floor(x/txw)	SSE4.1
						__m128 f_bch=_mm_floor_ps(bc_h);

						f_acw=_mm_mul_ps(f_acw, s_txwk);//floor(x/txw)*txw
						f_bch=_mm_mul_ps(f_bch, s_txhk);

						a_c=_mm_sub_ps(a_c, f_acw);//x-floor(x/txw)*txw
						b_c=_mm_sub_ps(b_c, f_bch);
						
						a_c=_mm_floor_ps(a_c);//99.9 -> 99 not 100	SSE4.1
						b_c=_mm_floor_ps(b_c);

						__m128i xtx=_mm_cvtps_epi32(a_c);
						__m128i ytx=_mm_cvtps_epi32(b_c);

						__m128i tx_idx=_mm_mullo_epi32(i_txwk, ytx);//SSE4.1
						tx_idx=_mm_add_epi32(tx_idx, xtx);

						__m128i m_tx=_mm_set_epi32(txk[tx_idx.m128i_i32[3]], txk[tx_idx.m128i_i32[2]], txk[tx_idx.m128i_i32[1]], txk[tx_idx.m128i_i32[0]]);
						__m128i m_rgb=_mm_loadu_si128((__m128i*)(rgb_y+xs));
						__m128i m_rgb2=_mm_and_si128(m_rgb, c_result2);
						m_tx=_mm_and_si128(m_tx, c_result);
						m_rgb=_mm_and_si128(m_rgb, c_result);

						const __m128i zero=_mm_setzero_si128();
						__m128i m_tx_lo=_mm_unpacklo_epi8(m_tx, zero);//{0, tx[7], ..., 0, tx[0]}
						__m128i m_tx_hi=_mm_unpackhi_epi8(m_tx, zero);//{0, tx[15], ..., 0, tx[8]}
						__m128i m_rgb_lo=_mm_unpacklo_epi8(m_rgb, zero);
						__m128i m_rgb_hi=_mm_unpackhi_epi8(m_rgb, zero);
						if(lsw_transparency_multiply)
						{
							m_tx_lo=_mm_mullo_epi16(m_tx_lo, m_rgb_lo);//{{0, tx[7]}*{0, rgb[7]}, ..., {0, tx[0]}*{0, rgb[0]}}
							m_tx_hi=_mm_mullo_epi16(m_tx_hi, m_rgb_hi);
							m_tx_lo=_mm_shuffle_epi8(m_tx_lo, sh_mul_lo);
							m_tx_hi=_mm_shuffle_epi8(m_tx_hi, sh_mul_hi);
						}
						else
						{
							m_tx_lo=_mm_add_epi16(m_tx_lo, m_rgb_lo);
							m_tx_hi=_mm_add_epi16(m_tx_hi, m_rgb_hi);
							m_tx_lo=_mm_srli_epi16(m_tx_lo, 1);
							m_tx_hi=_mm_srli_epi16(m_tx_hi, 1);
							m_tx_lo=_mm_shuffle_epi8(m_tx_lo, sh_add_lo);
							m_tx_hi=_mm_shuffle_epi8(m_tx_hi, sh_add_hi);
						}
						m_tx_lo=_mm_and_si128(m_tx_lo, mask_lo);
						m_tx_hi=_mm_and_si128(m_tx_hi, mask_hi);
						m_tx=_mm_or_si128(m_tx_lo, m_tx_hi);

						m_tx=_mm_or_si128(m_tx, m_rgb2);
						_mm_storeu_si128((__m128i*)(rgb_y+xs), m_tx);
					}
					adm=_mm_add_ps(adm, adm_increment);
					m_a=_mm_add_ps(m_a, m_a_increment);
					m_b=_mm_add_ps(m_b, m_b_increment);
					m_c=_mm_add_ps(m_c, m_c_increment);
				}
				admittance=adm.m128_f32[0];
				a=m_a.m128_f32[0];
				b=m_b.m128_f32[0];
				c=m_c.m128_f32[0];
				for(;xs<lfb_y;++xs)
				{
					if(std::abs(admittance-wb_y[xs])<1e-10*std::abs(admittance+wb_y[xs])||//z fighting
						admittance>wb_y[xs])
					{
						int Xtx=int(a/c)%txwk;Xtx+=txwk&-(Xtx<0);//if(Xtx<0)Xtx+=txwk;
						int Ytx=int(b/c)%txhk;Ytx+=txhk&-(Ytx<0);//if(Ytx<0)Ytx+=txhk;
						auto p=(unsigned char*)&rgb_y[xs], c=(unsigned char*)&txk[Xtx+Ytx*txwk];
						if(lsw_transparency_multiply)
						{
							p[0]=p[0]*c[0]>>8;//b
							p[1]=p[1]*c[1]>>8;//g
							p[2]=p[2]*c[2]>>8;//r
						}
						else
						{
							p[0]=(p[0]+c[0])>>1;//b
							p[1]=(p[1]+c[1])>>1;//g
							p[2]=(p[2]+c[2])>>1;//r
						}
					}
					admittance+=A1, a+=B1, b+=B4, c+=B7;
				}
				wb_y=wb_y+w, rgb_y=rgb_y+w;
			}
		}
		else//textured transparent triangle - no SSE4.1
		{
			const __m128 half=_mm_set1_ps(0.5f);
			const __m128i mask2_lo=_mm_set_epi32(0, -1, 0, -1), mask2_hi=_mm_set_epi32(-1, 0, -1, 0);
			const __m128i sh_mul_lo=_mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 15, 13, 11, 9, 7, 5, 3, 1);
			const __m128i sh_mul_hi=_mm_set_epi8(15, 13, 11, 9, 7, 5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0);
			const __m128i sh_add_lo=_mm_set_epi8(1, 1, 1, 1, 1, 1, 1, 1, 14, 12, 10, 8, 6, 4, 2, 0);
			const __m128i sh_add_hi=_mm_set_epi8(14, 12, 10, 8, 6, 4, 2, 0, 1, 1, 1, 1, 1, 1, 1, 1);
			const __m128i mask_lo=_mm_set_epi32(0, 0, -1, -1);
			const __m128i mask_hi=_mm_set_epi32(-1, -1, 0, 0);

			__m128 adm_increment=_mm_set1_ps(A1*4);
			__m128 m_a_increment=_mm_set1_ps(B1*4);
			__m128 m_b_increment=_mm_set1_ps(B4*4);
			__m128 m_c_increment=_mm_set1_ps(B7*4);

			__m128 s_txwk=_mm_set1_ps((float)txwk);
			__m128 s_txhk=_mm_set1_ps((float)txhk);
			__m128i i_txwk=_mm_set1_epi32(txwk), i_txwk_1=_mm_set1_epi32(txwk-1);
			__m128i i_txhk=_mm_set1_epi32(txhk), i_txhk_1=_mm_set1_epi32(txhk-1);
			for(int ys=0;ys<h;++ys)
			{
				lib_y=libuffer[ys]<0?0:libuffer[ys], lfb_y=lfbuffer[ys]>w?w:lfbuffer[ys];
				float admittance=A1*lib_y+A2*ys+A3, a=B1*lib_y+B2*ys+B3, b=B4*lib_y+B5*ys+B6, c=B7*lib_y+B8*ys+B9;

				int xs=lib_y;
				__m128 adm=_mm_set_ps(admittance+3*A1, admittance+A1+A1, admittance+A1, admittance);
				__m128 m_a=_mm_set_ps(a+B1*3, a+B1+B1, a+B1, a);
				__m128 m_b=_mm_set_ps(b+B4*3, b+B4+B4, b+B4, b);
				__m128 m_c=_mm_set_ps(c+B7*3, c+B7+B7, c+B7, c);
				for(int xsEnd=lfb_y-4;xs<xsEnd;xs+=4)//for each 4 pixels
				{
					__m128 wb=_mm_loadu_ps(wb_y+xs);
					__m128i c_result=_mm_castps_si128(_mm_cmplt_ps(wb, adm));
					__m128i c_result2=_mm_shuffle_epi32(c_result, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
					__m128i c_or=_mm_or_si128(c_result, c_result2);//r3|r1, r2|r0, r1|r3, r0|r2
					c_result2=_mm_shuffle_epi32(c_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
								c_or=_mm_or_si128(c_or, c_result2);//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
					if(c_or.m128i_i32[0])
					{
						c_result2=_mm_xor_si128(c_result, _mm_set1_epi32(-1));
						__m128 a_c=_mm_div_ps(m_a, m_c);
						__m128 b_c=_mm_div_ps(m_b, m_c);

						//mod operation		x%q = x-q*floor(x/q)
						__m128 ac_w=_mm_div_ps(a_c, s_txwk);
						__m128 bc_h=_mm_div_ps(b_c, s_txhk);

						//__m128 f_acw=_mm_floor_ps(ac_w);//floor(x/txw)	SSE4.1
						//__m128 f_bch=_mm_floor_ps(bc_h);
						ac_w=_mm_sub_ps(ac_w, half);//floor(x/q)
						bc_h=_mm_sub_ps(bc_h, half);
						__m128i i_acw=_mm_cvtps_epi32(ac_w);
						__m128i i_bch=_mm_cvtps_epi32(bc_h);
						__m128 f_acw=_mm_cvtepi32_ps(i_acw);
						__m128 f_bch=_mm_cvtepi32_ps(i_bch);

						f_acw=_mm_mul_ps(f_acw, s_txwk);//floor(x/txw)*txw
						f_bch=_mm_mul_ps(f_bch, s_txhk);

						a_c=_mm_sub_ps(a_c, f_acw);//x-floor(x/txw)*txw
						b_c=_mm_sub_ps(b_c, f_bch);
						
						//a_c=_mm_floor_ps(a_c);//99.9 -> 99 not 100	SSE4.1
						//b_c=_mm_floor_ps(b_c);
						//__m128i xtx=_mm_cvtps_epi32(a_c);
						//__m128i ytx=_mm_cvtps_epi32(b_c);
						a_c=_mm_sub_ps(a_c, half);
						b_c=_mm_sub_ps(b_c, half);
						__m128i xtx=_mm_cvtps_epi32(a_c);
						__m128i ytx=_mm_cvtps_epi32(b_c);

						__m128i cmp_x=_mm_cmplt_epi32(xtx, _mm_setzero_si128());//index guard
						__m128i cmp_y=_mm_cmplt_epi32(ytx, _mm_setzero_si128());
						cmp_x=_mm_and_si128(cmp_x, i_txwk);
						cmp_y=_mm_and_si128(cmp_y, i_txhk);
						xtx=_mm_add_epi32(xtx, cmp_x);
						ytx=_mm_add_epi32(ytx, cmp_y);
						cmp_x=_mm_cmpgt_epi32(xtx, i_txwk_1);
						cmp_y=_mm_cmpgt_epi32(ytx, i_txhk_1);
						cmp_x=_mm_and_si128(cmp_x, i_txwk);
						cmp_y=_mm_and_si128(cmp_y, i_txhk);
						xtx=_mm_sub_epi32(xtx, cmp_x);
						ytx=_mm_sub_epi32(ytx, cmp_y);

						//__m128i tx_idx=_mm_mullo_epi32(i_txwk, ytx);//SSE4.1
						__m128i tx_idx_lo=_mm_mul_epu32(ytx, i_txwk);
						ytx=_mm_shuffle_epi32(ytx, _MM_SHUFFLE(2, 3, 0, 1));//y2 y3 y0 y1
						__m128i tx_idx_hi=_mm_mul_epu32(ytx, i_txwk);
						tx_idx_hi=_mm_shuffle_epi32(tx_idx_hi, _MM_SHUFFLE(2, 3, 0, 1));
						tx_idx_lo=_mm_and_si128(tx_idx_lo, mask2_lo);
						tx_idx_hi=_mm_and_si128(tx_idx_hi, mask2_hi);
						__m128i tx_idx=_mm_or_si128(tx_idx_lo, tx_idx_hi);

						tx_idx=_mm_add_epi32(tx_idx, xtx);

						__m128i m_tx=_mm_set_epi32(txk[tx_idx.m128i_i32[3]], txk[tx_idx.m128i_i32[2]], txk[tx_idx.m128i_i32[1]], txk[tx_idx.m128i_i32[0]]);
						__m128i m_rgb=_mm_loadu_si128((__m128i*)(rgb_y+xs));
						__m128i m_rgb2=_mm_and_si128(m_rgb, c_result2);
						m_tx=_mm_and_si128(m_tx, c_result);
						m_rgb=_mm_and_si128(m_rgb, c_result);
						
						__m128i const zero=_mm_setzero_si128();
						__m128i m_tx_lo=_mm_unpacklo_epi8(m_tx, zero);//{0, tx[7], ..., 0, tx[0]}
						__m128i m_tx_hi=_mm_unpackhi_epi8(m_tx, zero);//{0, tx[15], ..., 0, tx[8]}
						__m128i m_rgb_lo=_mm_unpacklo_epi8(m_rgb, zero);
						__m128i m_rgb_hi=_mm_unpackhi_epi8(m_rgb, zero);
						if(lsw_transparency_multiply)
						{
							m_tx_lo=_mm_mullo_epi16(m_tx_lo, m_rgb_lo);//{{0, tx[7]}*{0, rgb[7]}, ..., {0, tx[0]}*{0, rgb[0]}}
							m_tx_hi=_mm_mullo_epi16(m_tx_hi, m_rgb_hi);
							m_tx_lo=_mm_shuffle_epi8(m_tx_lo, sh_mul_lo);
							m_tx_hi=_mm_shuffle_epi8(m_tx_hi, sh_mul_hi);
							m_tx_lo=_mm_and_si128(m_tx_lo, mask_lo);
							m_tx_hi=_mm_and_si128(m_tx_hi, mask_hi);
						}
						else
						{
							m_tx_lo=_mm_add_epi16(m_tx_lo, m_rgb_lo);
							m_tx_hi=_mm_add_epi16(m_tx_hi, m_rgb_hi);
							m_tx_lo=_mm_srli_epi16(m_tx_lo, 1);
							m_tx_hi=_mm_srli_epi16(m_tx_hi, 1);
							m_tx_lo=_mm_shuffle_epi8(m_tx_lo, sh_add_lo);
							m_tx_hi=_mm_shuffle_epi8(m_tx_hi, sh_add_hi);
						}
						m_tx=_mm_or_si128(m_tx_lo, m_tx_hi);

						m_tx=_mm_or_si128(m_tx, m_rgb2);
						_mm_storeu_si128((__m128i*)(rgb_y+xs), m_tx);
					}
					adm=_mm_add_ps(adm, adm_increment);
					m_a=_mm_add_ps(m_a, m_a_increment);
					m_b=_mm_add_ps(m_b, m_b_increment);
					m_c=_mm_add_ps(m_c, m_c_increment);
				}
				admittance=adm.m128_f32[0];
				a=m_a.m128_f32[0];
				b=m_b.m128_f32[0];
				c=m_c.m128_f32[0];
				for(;xs<lfb_y;++xs)
				{
					if(std::abs(admittance-wb_y[xs])<1e-10*std::abs(admittance+wb_y[xs])||//z fighting
						admittance>wb_y[xs])
					{
						int Xtx=int(a/c)%txwk;Xtx+=txwk&-(Xtx<0);//if(Xtx<0)Xtx+=txwk;
						int Ytx=int(b/c)%txhk;Ytx+=txhk&-(Ytx<0);//if(Ytx<0)Ytx+=txhk;
						auto p=(unsigned char*)&rgb_y[xs], c=(unsigned char*)&txk[Xtx+Ytx*txwk];
						if(lsw_transparency_multiply)
						{
							p[0]=p[0]*c[0]>>8;//b
							p[1]=p[1]*c[1]>>8;//g
							p[2]=p[2]*c[2]>>8;//r
						}
						else
						{
							p[0]=(p[0]+c[0])>>1;//b
							p[1]=(p[1]+c[1])>>1;//g
							p[2]=(p[2]+c[2])>>1;//r
						}
					}
					admittance+=A1, a+=B1, b+=B4, c+=B7;
				}
				wb_y=wb_y+w, rgb_y=rgb_y+w;
			}
		/*	for(int ys=0;ys<h;++ys)
			{
				lib_y=libuffer[ys]<0?0:libuffer[ys], lfb_y=lfbuffer[ys]>w?w:lfbuffer[ys];
				float admittance=A1*lib_y+A2*ys+A3, a=B1*lib_y+B2*ys+B3, b=B4*lib_y+B5*ys+B6, c=B7*lib_y+B8*ys+B9;
				for(int xs=lib_y;xs<lfb_y;++xs)
				{
					if(std::abs(admittance-wb_y[xs])<1e-10*std::abs(admittance+wb_y[xs])||//z fighting
						admittance>wb_y[xs])
					{
						int Xtx=int(a/c)%txwk;Xtx+=txwk&-(Xtx<0);//if(Xtx<0)Xtx+=txwk;
						int Ytx=int(b/c)%txhk;Ytx+=txhk&-(Ytx<0);//if(Ytx<0)Ytx+=txhk;
						auto p=(unsigned char*)&rgb_y[xs], c=(unsigned char*)&txk[Xtx+Ytx*txwk];
						p[0]=p[0]*c[0]>>8;//b
						p[1]=p[1]*c[1]>>8;//g
						p[2]=p[2]*c[2]>>8;//r
					}
					admittance+=A1, a+=B1, b+=B4, c+=B7;
				}
				wb_y=wb_y+w, rgb_y=rgb_y+w;
			}//*/
		}
	}
	else//solid color transparent triangle
	{
		float *wb_y=wbuffer;
		int lib_y, lfb_y, *rgb_y=rgb, color=tc;
		__m128i const sh_mul_lo=_mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 15, 13, 11, 9, 7, 5, 3, 1);
		__m128i const sh_mul_hi=_mm_set_epi8(15, 13, 11, 9, 7, 5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0);
		const __m128i sh_add_lo=_mm_set_epi8(1, 1, 1, 1, 1, 1, 1, 1, 14, 12, 10, 8, 6, 4, 2, 0);
		const __m128i sh_add_hi=_mm_set_epi8(14, 12, 10, 8, 6, 4, 2, 0, 1, 1, 1, 1, 1, 1, 1, 1);
		__m128i const mask_lo=_mm_set_epi32(0, 0, -1, -1);
		__m128i const mask_hi=_mm_set_epi32(-1, -1, 0, 0);
		__m128i m_color=_mm_set1_epi32(color);
		for(int ys=0;ys<h;++ys)
		{
			lib_y=libuffer[ys]<0?0:libuffer[ys], lfb_y=lfbuffer[ys]>w?w:lfbuffer[ys];
			float admittance=A1*lib_y+A2*ys+A3;

			int xs=lib_y;
			__m128 adm=_mm_set_ps(admittance+3*A1, admittance+A1+A1, admittance+A1, admittance);
			__m128 adm_increment=_mm_set1_ps(A1*4);
			for(int xsEnd=lfb_y-4;xs<xsEnd;xs+=4)
			{
				__m128 wb=_mm_loadu_ps(wb_y+xs);
				__m128i c_result=_mm_castps_si128(_mm_cmplt_ps(wb, adm));
				__m128i c_result2=_mm_shuffle_epi32(c_result, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
				__m128i c_or=_mm_or_si128(c_result, c_result2);//r3|r1, r2|r0, r1|r3, r0|r2
				c_result2=_mm_shuffle_epi32(c_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
				c_or=_mm_or_si128(c_or, c_result2);//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
				if(c_or.m128i_i32[0])
				{
					c_result2=_mm_xor_si128(c_result, _mm_set1_epi32(-1));

					__m128i m_rgb=_mm_loadu_si128((__m128i*)(rgb_y+xs));
					__m128i m_rgb2=_mm_and_si128(m_rgb, c_result2);
					__m128i m_c2=_mm_and_si128(m_color, c_result);
					m_rgb=_mm_and_si128(m_rgb, c_result);
						
					__m128i zero=_mm_setzero_si128();
					__m128i m_c_lo=_mm_unpacklo_epi8(m_c2, zero);//{0, tx[7], ..., 0, tx[0]}
					__m128i m_c_hi=_mm_unpackhi_epi8(m_c2, zero);//{0, tx[15], ..., 0, tx[8]}
					__m128i m_rgb_lo=_mm_unpacklo_epi8(m_rgb, zero);
					__m128i m_rgb_hi=_mm_unpackhi_epi8(m_rgb, zero);
					if(lsw_transparency_multiply)
					{
						m_c_lo=_mm_mullo_epi16(m_c_lo, m_rgb_lo);//{{0, tx[7]}*{0, rgb[7]}, ..., {0, tx[0]}*{0, rgb[0]}}
						m_c_hi=_mm_mullo_epi16(m_c_hi, m_rgb_hi);
						m_c_lo=_mm_shuffle_epi8(m_c_lo, sh_mul_lo);
						m_c_hi=_mm_shuffle_epi8(m_c_hi, sh_mul_hi);
					}
					else
					{
						m_c_lo=_mm_add_epi16(m_c_lo, m_rgb_lo);
						m_c_hi=_mm_add_epi16(m_c_hi, m_rgb_hi);
						m_c_lo=_mm_srli_epi16(m_c_lo, 1);
						m_c_hi=_mm_srli_epi16(m_c_hi, 1);
						m_c_lo=_mm_shuffle_epi8(m_c_lo, sh_add_lo);
						m_c_hi=_mm_shuffle_epi8(m_c_hi, sh_add_hi);
					}
					m_c_lo=_mm_and_si128(m_c_lo, mask_lo);
					m_c_hi=_mm_and_si128(m_c_hi, mask_hi);
					m_c2=_mm_or_si128(m_c_lo, m_c_hi);

					m_c2=_mm_or_si128(m_c2, m_rgb2);
					_mm_storeu_si128((__m128i*)(rgb_y+xs), m_c2);
				}
				adm=_mm_add_ps(adm, adm_increment);
			}
			admittance=adm.m128_f32[0];
			for(;xs<lfb_y;++xs)
			{
				if(std::abs(admittance-wb_y[xs])<1e-10*std::abs(admittance+wb_y[xs])||//z fighting
					admittance>=wb_y[xs])
				{
					auto p=(unsigned char*)&rgb_y[xs], c=(unsigned char*)&color;
					if(lsw_transparency_multiply)
					{
						p[0]=p[0]*c[0]>>8;//b
						p[1]=p[1]*c[1]>>8;//g
						p[2]=p[2]*c[2]>>8;//r
					}
					else
					{
						p[0]=(p[0]+c[0])>>1;//b
						p[1]=(p[1]+c[1])>>1;//g
						p[2]=(p[2]+c[2])>>1;//r
					}
				}
				admittance+=A1;
			}
			wb_y=wb_y+w, rgb_y=rgb_y+w;
		}
	}
//	GUIPrint(ghMemDC, (int)v1[0], (int)v1[1], k);
}
void				LinearSW::render()
{
#if 1
	for(int k=0;k<npols;++k)//draw opaque (ordinary) polygons first
	{
		auto &p=env[k];
		render_opaque(p.p1, p.p2, p.p3, p.tx_idx>=0, p.tx_idx>=0?p.tx_idx:p.color, p.tx1, p.txm);
	}
	for(int k=0, kEnd=objects.size();k<kEnd;++k)
	{
		auto &A=objects[k];
		for(int kt=0, ktEnd=A.tr.size();kt<ktEnd;++kt)
		{
			auto &tr=A.tr[kt];
			render_transparent(A.w[tr.a], A.w[tr.b], A.w[tr.c], tr.textured, tr.tc, tr.tx1, tr.txm);
		//	render_opaque(A.w[tr.a], A.w[tr.b], A.w[tr.c], tr.textured, tr.tc, tr.tx1, tr.txm);
		}
	}
	for(int k=0;k<ntpols;++k)//draw transparent triangles afterwards
	{
		auto &p=tenv[k];
		render_transparent(p.p1, p.p2, p.p3, p.tx_idx>=0, p.tx_idx>=0?p.tx_idx:p.color, p.tx1, p.txm);
	}//*/
	//Vector3f zero(0, 0, 0);
	//draw_line(zero, Vector3f(infinity, 0, 0), 0);
	//draw_line(zero, Vector3f(0, infinity, 0), 0);
	//draw_line(zero, Vector3f(0, 0, infinity), 0);
#endif

	QueryPerformanceFrequency(&li);freq=li.QuadPart;
	QueryPerformanceCounter(&li);
#if 1
		static float total_t, min_t=infinity;
		static int n_t=1900;//starting point
		if(++n_t>2000)
	//	if(++n_t>700)//less coherent?
	//	if(++n_t>100)
			n_t=1, total_t=0, min_t=infinity;
		float dt=1000.f*(li.QuadPart-nticks)/freq;
		if(min_t>dt)
			min_t=dt;
		total_t+=dt;
		SetBkMode(ghMemDC, OPAQUE);
		GUIPrint(ghMemDC, 0, 0, "T_min=%.10f, T_av=%.10f, NT=%d", min_t, total_t/n_t, n_t);
		SetBkMode(ghMemDC, TRANSPARENT);
#endif
	if(w&&h)
		linelen=sprintf_s(line, 128, "fps=%.10f, T=%.10fms, dcam=%.10f, fov=%.10f, distance=%g", freq/float(li.QuadPart-nticks), 1000.*(li.QuadPart-nticks)/freq, dcam, atan(tanfov)*720/_2pi, 1/wbuffer[w*(h/2)+w/2]);
	else
		linelen=sprintf_s(line, 128, "fps=%.10f, T=%.10fms, dcam=%.10f, fov=%.10f", freq/float(li.QuadPart-nticks), 1000.*(li.QuadPart-nticks)/freq, dcam, atan(tanfov)*720/_2pi);
	nticks=li.QuadPart;
	if(linelen>0)TextOutA(ghMemDC, 0, h-16, line, linelen);
//	GUIPrint(ghMemDC, 0,  0, "lParam=0x%08X", glParam);
//	GUIPrint(ghMemDC, 0, 18, "wParam=0x%08X", gwParam);
}
void				LinearSW::show(){BitBlt(ghDC, 0, 0, w, h, ghMemDC, 0, 0, SRCCOPY);}
void				LinearSW::draft_1behind		(Vector2f &s1, Vector2f &s2, Vector2f &s3)//1 & 2 forward		3 behind
{
	draft_start(s1, s2);
	draft_crit(s3, s1);
	draft_crit(s3, s2);
	if(s3.y>s1.y&&s3.y<s2.y||s3.y>s2.y&&s3.y<s1.y){		 if(s3.x<(s2.x-s1.x)*(s3.y-s1.y)/(s2.y-s1.y)+s1.x)	for(int k=0;k<h;++k)lfbuffer[k]=w;
													else													memset(libuffer, 0, h*sizeof(int));}
}
void				LinearSW::draft_2behind		(Vector2f &s1, Vector2f &s2, Vector2f &s3)//1 forward		2 & 3 behind
{
	draft_crit_start(s2, s1);
	draft_crit(s3, s1);
	if(s1.y>s2.y&&s1.y<s3.y||s1.y>s3.y&&s1.y<s2.y){		 if(s1.x<(s3.x-s2.x)*(s1.y-s2.y)/(s3.y-s2.y)+s2.x)	memset(libuffer, 0, h*sizeof(int));
													else													for(int k=0;k<h;++k)lfbuffer[k]=w;}
}
void				LinearSW::draft_start		(Vector2f &s1, Vector2f &s2)
{
	for(int k2=0;k2<h;++k2)libuffer[k2]=w; memset(lfbuffer, 0, h*sizeof(int));
	int k2;
	float k3, A=(s2.x-s1.x)/(s2.y-s1.y);
	if(s1.y<s2.y)
	{
		k3=s1.x+A*((long long)clamp_positive(s1.y)+1-s1.y);
		//k3=s1.x+A*((long long)(s1.y<0?0:s1.y)+1-s1.y);
		if(s1.x<s2.x)
			for(long long k=long long(s1.y)<0?	0:long long(s1.y)	+1;k<h&&k<s2.y	;++k)
				 k2=int(k3), k2=k2<0
					?	(int)clamp_positive(s1.x)
					:	k2>w
						?	minimum(int(s2.x), w)
						:	k2<s1.x
							?	int(s1.x)
							:	minimum(int(s2.x), k2), libuffer[k]=lfbuffer[k]=k2, k3+=A;
			//	 k2=int(k3), k2=k2<0?s1.x>0?int(s1.x):0:k2>w?s2.x<w?int(s2.x):w:k2<s1.x?int(s1.x):k2>s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;
		else
			for(long long k=long long(s1.y)<0?	0:long long(s1.y)	+1;k<h&&k<s2.y	;++k)
				k2=int(k3), k2=k2<0
					?	(int)clamp_positive(s2.x)
					:	k2>w
						?	minimum(int(s1.x), w)
						:	k2<s2.x
							?	int(s2.x)
							:	minimum(int(s1.x), k2), libuffer[k]=lfbuffer[k]=k2, k3+=A;
			//	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?s1.x<w?int(s1.x):w:k2<s2.x?int(s2.x):k2>s1.x?int(s1.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;
	}
	else
	{
		k3=s1.x+A*((long long)(s2.y<0?0:s2.y)+1-s1.y);
		if(s1.x<s2.x)
		{
			int kEnd=minimum(int(ceil(s1.y)), h);
			for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<kEnd	;++k)
		//	for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h&&k<s1.y	;++k)
				k2=int(k3), k2=k2<0
					?	(int)clamp_positive(s1.x)
					:	k2>w
						?	minimum(int(s2.x), w)
						:	k2<s1.x
							?	int(s1.x)
							:	minimum(int(s2.x), k2), libuffer[k]=lfbuffer[k]=k2, k3+=A;
			//	k2=int(k3), k2=k2<0?s1.x>0?int(s1.x):0:k2>w?s2.x<w?int(s2.x):w:k2<s1.x?int(s1.x):k2>s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;
		}
		else
		{
			int kEnd=minimum(int(ceil(s1.y)), h);
			for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<kEnd	;++k)
		//	for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h&&k<s1.y	;++k)
				k2=int(k3), k2=k2<0
					?	(int)clamp_positive(int(s2.x))
					:	k2>w
						?	minimum(int(s1.x), w)
						:	k2<s2.x
							?	int(s2.x)
							:	minimum(int(s1.x), k2), libuffer[k]=lfbuffer[k]=k2, k3+=A;
			//	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?s1.x<w?int(s1.x):w:k2<s2.x?int(s2.x):k2>s1.x?int(s1.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;
		}
	}
	//	 if(s1.y<s2.y){	k3=s1.x+A*((long long)(s1.y<0?0:s1.y)+1-s1.y);			 if(s1.x<s2.x)	for(long long k=long long(s1.y)<0?	0:long long(s1.y)	+1;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?s1.x>0?int(s1.x):0:k2>w?s2.x<w?int(s2.x):w:k2<s1.x?int(s1.x):k2>s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;
	//																		else				for(long long k=long long(s1.y)<0?	0:long long(s1.y)	+1;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?s1.x<w?int(s1.x):w:k2<s2.x?int(s2.x):k2>s1.x?int(s1.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;}
	//else				{	k3=s1.x+A*((long long)(s2.y<0?0:s2.y)+1-s1.y);		 if(s1.x<s2.x)	for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h&&k<s1.y	;++k)	k2=int(k3), k2=k2<0?s1.x>0?int(s1.x):0:k2>w?s2.x<w?int(s2.x):w:k2<s1.x?int(s1.x):k2>s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;
	//																		else				for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h&&k<s1.y	;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?s1.x<w?int(s1.x):w:k2<s2.x?int(s2.x):k2>s1.x?int(s1.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;}
}
void				LinearSW::draft				(Vector2f &s1, Vector2f &s2)
{
	int k2;float k3, A=(s2.x-s1.x)/(s2.y-s1.y);
	if(s1.y<s2.y)
	{
		k3=s1.x+A*((long long)(s1.y<0?0:s1.y)+1-s1.y);
		if(s1.x<s2.x)
		{
			int kEnd=minimum(int(ceil(s2.y)), h);
			for(long long k=long long(s1.y)<0?	0:long long(s1.y)	+1;k<kEnd	;++k)
		//	for(long long k=long long(s1.y)<0?	0:long long(s1.y)	+1;k<h&&k<s2.y	;++k)
				k2=int(k3), k2=k2<0
					?	(int)clamp_positive(s1.x)
					:	k2>w
						?	minimum(int(s2.x), w)
						:	k2<s1.x
							?	int(s1.x)
							:	minimum(int(s2.x), k2), libuffer[k]=minimum(libuffer[k], k2), lfbuffer[k]=maximum(lfbuffer[k], k2), k3+=A;
		}
		else
		{
			int kEnd=minimum(int(ceil(s2.y)), h);
			for(long long k=long long(s1.y)<0?	0:long long(s1.y)	+1;k<kEnd	;++k)
		//	for(long long k=long long(s1.y)<0?	0:long long(s1.y)	+1;k<h&&k<s2.y	;++k)
				k2=int(k3), k2=k2<0
					?	(int)clamp_positive(s2.x)
					:	k2>w
						?	minimum(int(s1.x), w)
						:	k2<s2.x
							?	int(s2.x)
							:	minimum(int(s1.x), k2), libuffer[k]=minimum(libuffer[k], k2), lfbuffer[k]=maximum(lfbuffer[k], k2), k3+=A;
		}
	}
	else
	{
		k3=s1.x+A*((long long)(s2.y<0?0:s2.y)+1-s1.y);
		if(s1.x<s2.x)
		{
			int kEnd=minimum(int(ceil(s1.y)), h);
			for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<kEnd	;++k)
		//	for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h&&k<s1.y	;++k)
				k2=int(k3), k2=k2<0
					?	(int)clamp_positive(s1.x)
					:	k2>w
						?	minimum(int(s2.x), w)
						:	k2<s1.x
							?	int(s1.x)
							:	minimum(int(s2.x), k2), libuffer[k]=minimum(libuffer[k], k2), lfbuffer[k]=maximum(lfbuffer[k], k2), k3+=A;
		}
		else
		{
			int kEnd=minimum(int(ceil(s1.y)), h);
			for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<kEnd	;++k)
		//	for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h&&k<s1.y	;++k)
				k2=int(k3), k2=k2<0
					?	(int)clamp_positive(s2.x)
					:	k2>w
						?	minimum(int(s1.x), w)
						:	k2<s2.x
							?	int(s2.x)
							:	minimum(int(s1.x), k2), libuffer[k]=minimum(libuffer[k], k2), lfbuffer[k]=maximum(lfbuffer[k], k2), k3+=A;
		}
	}
	//	 if(s1.y<s2.y)	{k3=s1.x+A*((long long)(s1.y<0?0:s1.y)+1-s1.y);			 if(s1.x<s2.x)	for(long long k=long long(s1.y)<0?	0:long long(s1.y)	+1;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?s1.x>0?int(s1.x):0:k2>w?s2.x<w?int(s2.x):w:k2<s1.x?int(s1.x):k2>s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;
	//																		else				for(long long k=long long(s1.y)<0?	0:long long(s1.y)	+1;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?s1.x<w?int(s1.x):w:k2<s2.x?int(s2.x):k2>s1.x?int(s1.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;}
	//else				{k3=s1.x+A*((long long)(s2.y<0?0:s2.y)+1-s1.y);			 if(s1.x<s2.x)	for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h&&k<s1.y	;++k)	k2=int(k3), k2=k2<0?s1.x>0?int(s1.x):0:k2>w?s2.x<w?int(s2.x):w:k2<s1.x?int(s1.x):k2>s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;
	//																		else				for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h&&k<s1.y	;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?s1.x<w?int(s1.x):w:k2<s2.x?int(s2.x):k2>s1.x?int(s1.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;}
}
void				LinearSW::draft_crit_start	(Vector2f &s1, Vector2f &s2)
{
	for(int k2=0;k2<h;++k2)libuffer[k2]=w; memset(lfbuffer, 0, h*sizeof(int));
	int k2;float k3, A=(s2.x-s1.x)/(s2.y-s1.y);
	if(s1.y<s2.y)
	{
		k3=s1.x+A*((long long)(s2.y<0?0:s2.y)+1-s1.y);
		if(s1.x<s2.x)
			for(long long k=long long(s2.y)<0?	0:long long(s2.y)+1	;k<h			;++k)
				k2=int(k3), k2=k2<0
					?	(int)clamp_positive(s2.x)
					:	k2>w
						?	w
						:	maximum(int(s2.x), k2), libuffer[k]=lfbuffer[k]=k2, k3+=A;
		else
			for(long long k=long long(s2.y)<0?	0:long long(s2.y)+1	;k<h			;++k)
				k2=int(k3), k2=k2<0
					?	0
					:	minimum(int(s2.x), k2, w), libuffer[k]=lfbuffer[k]=k2, k3+=A;
				//k2=int(k3), k2=k2<0
				//	?	0
				//	:	k2>w
				//		?	minimum(int(s2.x), w)
				//		:	minimum(int(s2.x), k2), libuffer[k]=lfbuffer[k]=k2, k3+=A;
	}
	else
	{
		k3=s1.x-A*((long long)s1.y+1);
		if(s1.x<s2.x)
			for(long long k=					0					, kEnd=minimum(int(ceil(s2.y)), h);k<kEnd	;++k)
		//	for(long long k=					0					;k<h&&k<s2.y	;++k)
				k2=int(k3), k2=k2<0
					?	(int)clamp_positive(s2.x)
					:	k2>w
						?	w
						:	maximum(int(s2.x), k2), libuffer[k]=lfbuffer[k]=k2, k3+=A;
		else
			for(long long k=					0					, kEnd=minimum(int(ceil(s2.y)), h);k<kEnd	;++k)
		//	for(long long k=					0					;k<h&&k<s2.y	;++k)
				k2=int(k3), k2=k2<0
					?	0
					:	minimum(int(s2.x), k2, w), libuffer[k]=lfbuffer[k]=k2, k3+=A;
	}
	//	 if(s1.y<s2.y)	{	k3=s1.x+A*((long long)(s2.y<0?0:s2.y)+1-s1.y);		 if(s1.x<s2.x)	for(long long k=long long(s2.y)<0?	0:long long(s2.y)+1	;k<h			;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?w:k2<s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;
	//																		else				for(long long k=long long(s2.y)<0?	0:long long(s2.y)+1	;k<h			;++k)	k2=int(k3), k2=k2<0?0:k2>w?s2.x<w?int(s2.x):w:k2>s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;}
	//else				{	k3=s1.x-A*((long long)s1.y+1);						 if(s1.x<s2.x)	for(long long k=					0					;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?w:k2<s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;
	//																		else				for(long long k=					0					;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?0:k2>w?s2.x<w?int(s2.x):w:k2>s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;}

	//	 if(s1.y<s2.y)	{	k3=s1.x+A*((long long)(s2.y<0?0:s2.y)+1-s1.y);		 if(s1.x<s2.x)	for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h			;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?w:k2<s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;
	//																		else				for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h			;++k)	k2=int(k3), k2=k2<0?0:k2>w?s2.x<w?int(s2.x):w:k2>s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;}
	//else				{	k3=s1.x-A*((long long)s1.y+1);						 if(s1.x<s2.x)	for(long long k=					0					+1;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?w:k2<s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;
	//																		else				for(long long k=					0					+1;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?0:k2>w?s2.x<w?int(s2.x):w:k2>s2.x?int(s2.x):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;}
}
void				LinearSW::draft_crit		(Vector2f &s1, Vector2f &s2)
{
	int k2;float k3, A=(s2.x-s1.x)/(s2.y-s1.y);
	if(s1.y<s2.y)
	{
		k3=s1.x+A*((long long)(s2.y<0?0:s2.y)+1-s1.y);
		if(s1.x<s2.x)
			for(long long k=long long(s2.y)<0?	0:long long(s2.y)+1	;k<h			;++k)
				k2=int(k3), k2=k2<0
					?	(int)clamp_positive(s2.x)
					:	k2>w
						?	w
						:	maximum(int(s2.x), k2), libuffer[k]=minimum(libuffer[k], k2), lfbuffer[k]=maximum(lfbuffer[k], k2), k3+=A;
		else
			for(long long k=long long(s2.y)<0?	0:long long(s2.y)+1	;k<h			;++k)
				k2=int(k3), k2=k2<0
					?	0
					:	minimum(int(s2.x), k2, w), libuffer[k]=minimum(libuffer[k], k2), lfbuffer[k]=maximum(lfbuffer[k], k2), k3+=A;
	}
	else
	{
		k3=s1.x-A*((long long)s1.y+1);
		if(s1.x<s2.x)
			for(long long k=					0					, kEnd=minimum(int(s2.y)+1, h);k<kEnd	;++k)
		//	for(long long k=					0					;k<h&&k<s2.y	;++k)
				k2=int(k3), k2=k2<0
					?	(int)clamp_positive(s2.x)
					:	k2>w
						?	w
						:	maximum(int(s2.x), k2), libuffer[k]=minimum(libuffer[k], k2), lfbuffer[k]=maximum(lfbuffer[k], k2), k3+=A;
		else
			for(long long k=					0					, kEnd=minimum(int(s2.y)+1, h);k<kEnd	;++k)
		//	for(long long k=					0					;k<h&&k<s2.y	;++k)
				k2=int(k3), k2=k2<0
					?	0
					:	minimum(int(s2.x), k2, w), libuffer[k]=minimum(libuffer[k], k2), lfbuffer[k]=maximum(lfbuffer[k], k2), k3+=A;
	}
	//	 if(s1.y<s2.y)	{	k3=s1.x+A*((long long)(s2.y<0?0:s2.y)+1-s1.y);		 if(s1.x<s2.x)	for(long long k=long long(s2.y)<0?	0:long long(s2.y)+1	;k<h			;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?w:k2<s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;
	//																		else				for(long long k=long long(s2.y)<0?	0:long long(s2.y)+1	;k<h			;++k)	k2=int(k3), k2=k2<0?0:k2>w?s2.x<w?int(s2.x):w:k2>s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;}
	//else				{	k3=s1.x-A*((long long)s1.y+1);						 if(s1.x<s2.x)	for(long long k=					0					;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?w:k2<s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;
	//																		else				for(long long k=					0					;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?0:k2>w?s2.x<w?int(s2.x):w:k2>s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;}
	//if(s1.y<s2.y)
	//{
	//	k3=s1.x+A*((long long)(s2.y<0?0:s2.y)+1-s1.y);
	//	if(s1.x<s2.x)
	//		for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h			;++k)
	//			k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?w:k2<s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;
	//	else
	//		for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h			;++k)
	//			k2=int(k3), k2=k2<0?0:k2>w?s2.x<w?int(s2.x):w:k2>s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;
	//}
	//else
	//{
	//	k3=s1.x-A*((long long)s1.y+1);
	//	if(s1.x<s2.x)for(long long k=			0					+1;k<h&&k<s2.y	;++k)
	//		k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?w:k2<s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;
	//	else
	//		for(long long k=					0					+1;k<h&&k<s2.y	;++k)
	//			k2=int(k3), k2=k2<0?0:k2>w?s2.x<w?int(s2.x):w:k2>s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;
	//}
//		 if(s1.y<s2.y){	k3=s1.x+A*((long long)(s2.y<0?0:s2.y)+1-s1.y);		 if(s1.x<s2.x)for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h			;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?w:k2<s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;
//																				else				for(long long k=long long(s2.y)<0?	0:long long(s2.y)	+1;k<h			;++k)	k2=int(k3), k2=k2<0?0:k2>w?s2.x<w?int(s2.x):w:k2>s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;}
//	else				{	k3=s1.x-A*((long long)s1.y+1);						 if(s1.x<s2.x)for(long long k=					0					+1;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?s2.x>0?int(s2.x):0:k2>w?w:k2<s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;
//																				else				for(long long k=					0					+1;k<h&&k<s2.y	;++k)	k2=int(k3), k2=k2<0?0:k2>w?s2.x<w?int(s2.x):w:k2>s2.x?int(s2.x):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;}
}

void				LinearSW::pushTexture(int *texture_)
{
	++ntx;
	texture=(int**)realloc(texture, ntx*sizeof(int*)), texture[ntx-1]=texture_;
}
void				LinearSW::enqueueTextures(){}
int*				LinearSW::popTexture()
{
	ntx--;
	return texture[ntx];
}
void				LinearSW::clearTextures()
{
	for(int k=0;k<ntx;++k)
		free(texture[k]);
	ntx=0;
}
void				LinearSW::finish()
{
	DeleteObject(SelectObject(ghMemDC, hBitmap)), DeleteDC(ghMemDC);
	for(int k=0;k<ntx;++k)free(texture[k]);
	free(texture), free(libuffer), free(lfbuffer), free(wbuffer);
//	free(rgb_r), free(rgb_g), free(rgb_b), free(rgb_n);

//	free(container);
}
struct			LinearGL:private Mode
{
	unsigned int *texture;
	HGLRC__ *hRC;
	void initiate();
	void resize();
	void changeFov();
	void clear_screen();
	void print(int x, int y, const char *a, ...);
	void print(Vector3f const &p, const char *a, ...);
	void draw_point(Vector3f const &p, int color);
	void draw_line(float x1, float y1, float x2, float y2, int lineColor);
	void draw_line(Vector3f const &p1, Vector3f const &p2, int lineColor);
	void draw_triangle(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, int color, bool transparent);
	void draw_ground();
	void render_opaque		(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, bool textured, int tc, Vector2f const &tx1, Vector2f const &tx2, Vector2f const &tx3);
	void render_transparent	(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, bool textured, int tc, Vector2f const &tx1, Vector2f const &tx2, Vector2f const &tx3);
	void render();
	void show();
	void pushTexture(int *texture_);
	void enqueueTextures();
	int* popTexture();
	void clearTextures();
	void finish();
} lgl;
void				LinearGL::initiate()
{
	texture=0;
	tagPIXELFORMATDESCRIPTOR pfd={sizeof(tagPIXELFORMATDESCRIPTOR), 1, PFD_DRAW_TO_WINDOW|PFD_SUPPORT_OPENGL|PFD_DOUBLEBUFFER, PFD_TYPE_RGBA, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, PFD_MAIN_PLANE, 0, 0, 0, 0};
	int PixelFormat=ChoosePixelFormat(ghDC, &pfd);
	SetPixelFormat(ghDC, PixelFormat, &pfd);
	hRC=wglCreateContext(ghDC);
	wglMakeCurrent(ghDC, hRC);
	if(h<=0)h=1;
	changeFov();
	glShadeModel(GL_SMOOTH);
	glClearColor(1, 1, 1, 0);
	glClearDepth(1);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	if(!glWindowPos2i)
		glWindowPos2i=(void(__stdcall*)(int, int))glutGetProcAddress("glWindowPos2i");
}
void				LinearGL::resize()
{
	if(h<=0)h=1;
	changeFov();
	SetBkMode(ghDC, OPAQUE);
}
void				LinearGL::changeFov()
{
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(atan(tanfov*h/w)*720/_2pi, float(w)/h, 50, 50000);//
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}
void				LinearGL::clear_screen(){glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);}
void				LinearGL::print(int x, int y, const char *a, ...)
{
	g_buflen=vsprintf_s(g_buf, g_bufsize, a, (char*)(&a+1));

	glColor3f(0, 0, 0);
	glWindowPos2i(x, h-18-y);
//	glRasterPos2f((float)x, (float)y);
	for(int k=0;k<g_buflen;++k)
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, g_buf[k]);
}
void				LinearGL::print(Vector3f const &p, const char *a, ...)
{//https://stackoverflow.com/questions/2183270/what-is-the-easiest-way-to-print-text-to-screen-in-opengl
	Vector2f s;
	Vector3f c;
	world_to_screen(p, s, c);
	if(c.z>0)
	{
		g_buflen=vsprintf_s(g_buf, g_bufsize, a, (char*)(&a+1));

		glColor3f(0, 0, 0);
		glWindowPos2i((int)s.x, h-(int)s.y);
	//	glRasterPos2f(s.x, s.y);
		for(int k=0;k<g_buflen;++k)
			glutBitmapCharacter(GLUT_BITMAP_9_BY_15, g_buf[k]);
	}
}
void				LinearGL::draw_point(Vector3f const &p, int color)
{
	Vector3f cp;
	world_to_camera(p, cp);

	glBegin(GL_POINTS);
	glVertex3f(cp.x, -cp.y, -cp.z);
	glEnd();
}
void				LinearGL::draw_line(float x1, float y1, float x2, float y2, int lineColor)
{
	//Vector3f zero(0, 0, 0), xy(1000, 1000, 0);
	//Vector3f zero_c, x_c;
	//world_to_camera(zero, zero_c), world_to_camera(xy, x_c);

	glBegin(GL_LINES);
	glColor3f(float(((unsigned char*)&lineColor)[2])/255, float(((unsigned char*)&lineColor)[1])/255, float(((unsigned char*)&lineColor)[0])/255);
	glVertex3f((x1-X0)*0.05f, (Y0-y1)*0.05f, -50);
	glVertex3f((x2-X0)*0.05f, (Y0-y2)*0.05f, -50);
	//glVertex3f(x1-X0, Y0-y1, -1000);
	//glVertex3f(x2-X0, Y0-y2, -1000);
	//glVertex3f(zero_c.x, -zero_c.y, -zero_c.z);
	//glVertex3f(x_c.x, -x_c.y, -x_c.z);
	//glVertex2f(x1, y1);
	//glVertex2f(x2, y2);
	glEnd();
}
void				LinearGL::draw_line(Vector3f const &p1, Vector3f const &p2, int lineColor)
{
	Vector3f c1, c2;
	world_to_camera(p1, c1), world_to_camera(p2, c2);

//	glDisable(GL_TEXTURE_2D);
	glBegin(GL_LINES);
	glColor3f(float(((unsigned char*)&lineColor)[2])/255, float(((unsigned char*)&lineColor)[1])/255, float(((unsigned char*)&lineColor)[0])/255);
	glVertex3f(c1.x, -c1.y, -c1.z);
	glVertex3f(c2.x, -c2.y, -c2.z);
	glEnd();
}
void				LinearGL::draw_triangle(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, int color, bool transparent)
{
	Vector2f zero;
	if(transparent)
		render_transparent(p1, p2, p3, false, color, zero, zero, zero);
	else
		render_opaque(p1, p2, p3, false, color, zero, zero, zero);
}
void				LinearGL::draw_ground()
{
}
void				LinearGL::render_opaque(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, bool textured, int tc, Vector2f const &tx1, Vector2f const &tx2, Vector2f const &tx3)
{
	Vector3f c1, c2, c3;
	world_to_camera(p1, c1), world_to_camera(p2, c2), world_to_camera(p3, c3);
//	float v1[3], v2[3], v3[3];//*v: Xcp, Ycp, Zcp
//	world_to_camera(p1, v1), world_to_camera(p2, v2), world_to_camera(p3, v3);
	if(textured)
	{
		glEnable(GL_TEXTURE_2D);
		int txwk=txw[tc], txhk=txh[tc];
		glBindTexture(GL_TEXTURE_2D, texture[tc]);
		glBegin(GL_TRIANGLES);//glDrawArrays
		glColor3f(1, 1, 1);
		glTexCoord2f(tx1.x/txwk, tx1.y/txhk), glVertex3f(c1.x, -c1.y, -c1.z);
		glTexCoord2f(tx2.x/txwk, tx2.y/txhk), glVertex3f(c2.x, -c2.y, -c2.z);
		glTexCoord2f(tx3.x/txwk, tx3.y/txhk), glVertex3f(c3.x, -c3.y, -c3.z);
		glEnd();
	}
	else
	{
		glDisable(GL_TEXTURE_2D);
		glBegin(GL_TRIANGLES);
		glColor3f(float(((unsigned char*)&tc)[2])/255, float(((unsigned char*)&tc)[1])/255, float(((unsigned char*)&tc)[0])/255);
	//	glColor3f(float(unsigned char(tc>>16))/255, float(unsigned char(tc>>8))/255, float(unsigned char(tc))/255);
		glVertex3f(c1.x, -c1.y, -c1.z);
		glVertex3f(c2.x, -c2.y, -c2.z);
		glVertex3f(c3.x, -c3.y, -c3.z);
		glEnd();
	}
}
void				LinearGL::render_transparent(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, bool textured, int tc, Vector2f const &tx1, Vector2f const &tx2, Vector2f const &tx3)
{
	Vector3f c1, c2, c3;
	world_to_camera(p1, c1), world_to_camera(p2, c2), world_to_camera(p3, c3);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	if(textured)
	{
		glEnable(GL_TEXTURE_2D);
		int txwk=txw[tc], txhk=txh[tc];
		glBindTexture(GL_TEXTURE_2D, texture[tc]);
		glBegin(GL_TRIANGLES);//glDrawArrays
		glColor4f(1, 1, 1, 0.5f);
	//	glColor3f(1, 1, 1);
		glTexCoord2f(tx1.x/txwk, tx1.y/txhk), glVertex3f(c1.x, -c1.y, -c1.z);
		glTexCoord2f(tx2.x/txwk, tx2.y/txhk), glVertex3f(c2.x, -c2.y, -c2.z);
		glTexCoord2f(tx3.x/txwk, tx3.y/txhk), glVertex3f(c3.x, -c3.y, -c3.z);
		glEnd();
	}
	else
	{
		glDisable(GL_TEXTURE_2D);
		glBegin(GL_TRIANGLES);
		glColor4f(float(((unsigned char*)&tc)[2])/255, float(((unsigned char*)&tc)[1])/255, float(((unsigned char*)&tc)[0])/255, 0.5f);
	//	glColor3f(float(unsigned char(tc>>16))/255, float(unsigned char(tc>>8))/255, float(unsigned char(tc))/255);
		glVertex3f(c1.x, -c1.y, -c1.z);
		glVertex3f(c2.x, -c2.y, -c2.z);
		glVertex3f(c3.x, -c3.y, -c3.z);
		glEnd();
	}
}
void 				LinearGL::render()
{
	glLoadIdentity();
	glTranslatef(0, 0, 0);
	for(int k=0, kEnd=objects.size();k<kEnd;++k)
	{
		auto &A=objects[k];
		for(int kt=0, ktEnd=A.tr.size();kt<ktEnd;++kt)
		{
			auto &tr=A.tr[kt];
			render_transparent(A.w[tr.a], A.w[tr.b], A.w[tr.c], tr.textured, tr.tc, tr.tx1, tr.tx2, tr.tx3);
		//	render_opaque(A.w[tr.a], A.w[tr.b], A.w[tr.c], tr.textured, tr.tc, tr.tx1, tr.tx2, tr.tx3);
		}
	}
	for(int k=0;k<npols;++k)
	{
		auto &p=env[k];
		render_opaque(p.p1, p.p2, p.p3, p.tx_idx>=0, p.tx_idx>=0?p.tx_idx:p.color, p.tx1, p.tx2, p.tx3);
	}
	for(int k=0;k<ntpols;++k)
	{
		auto &p=tenv[k];
		render_transparent(p.p1, p.p2, p.p3, p.tx_idx>=0, p.tx_idx>=0?p.tx_idx:p.color, p.tx1, p.tx2, p.tx3);
	}
	//Vector3f origin(0, 0, 0);
	//draw_line(origin, Vector3f(10000, 0, 0), 0);
	//draw_line(origin, Vector3f(0, 10000, 0), 0);
	//draw_line(origin, Vector3f(0, 0, 10000), 0);

	QueryPerformanceFrequency(&li);freq=li.QuadPart;
	QueryPerformanceCounter(&li);
	linelen=sprintf_s(line, "OpenGL, fps=%.10f, T=%.10fms, dcam=%.10f, fov=%.10f", freq/float(li.QuadPart-nticks), 1000.*(li.QuadPart-nticks)/freq, dcam, atan(tanfov)*720/_2pi);
	nticks=li.QuadPart;
}
void				LinearGL::show()
{
	SwapBuffers(ghDC);
	if(linelen>0)
		TextOutA(ghDC, 0, R.bottom-16, line, linelen);
}
void				LinearGL::pushTexture(int *texture_)
{
	++ntx;
	for(int k=0;k<txh[ntx-1]*txw[ntx-1];++k)texture_[k]=unsigned char(texture_[k])<<16|unsigned char(texture_[k]>>8)<<8|unsigned char(texture_[k]>>16);
	texture=(unsigned int*)realloc(texture, ntx*sizeof(unsigned int));
	glGenTextures(1, &texture[ntx-1]);
	glBindTexture(GL_TEXTURE_2D, texture[ntx-1]);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, txw[ntx-1], txh[ntx-1], 0, GL_RGBA, GL_UNSIGNED_BYTE, texture_);
	free(texture_);
}
void				LinearGL::enqueueTextures()
{
}
int*				LinearGL::popTexture()
{
	int *texture_=(int*)malloc(txh[ntx-1]*txw[ntx-1]*sizeof(int));
	glBindTexture(GL_TEXTURE_2D, texture[ntx-1]);
	glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE, texture_);
	glDeleteTextures(1, &texture[ntx-1]);
	for(int k=0;k<txh[ntx-1]*txw[ntx-1];++k)texture_[k]=unsigned char(texture_[k])<<16|unsigned char(texture_[k]>>8)<<8|unsigned char(texture_[k]>>16);
	ntx--;
	return texture_;
}
void				LinearGL::clearTextures()
{
	glDeleteTextures(ntx, texture);
	ntx=0;
}
void				LinearGL::finish()
{
	wglMakeCurrent(0, 0);
	wglDeleteContext(hRC);
	free(texture);
}
struct			NonlinearSW:private Mode
{
	int **texture;
	float *wbuffer, *xbuffer, *ybuffer, *zbuffer;
	int *rgb, rgbn;
	HDC__ *ghMemDC;
	HBITMAP__ *hBitmap;
	void (NonlinearSW::*rebuffer)();
	void initiate();
	void resize();
	void changeFov();
	void clear_screen();
	void print(int x, int y, const char *a, ...);
	void print(Vector3f const &p, const char *a, ...);
	void draw_point(Vector3f const &p, int color);
	void draw_line(float x1, float y1, float x2, float y2, int lineColor);
	void draw_line(Vector3f const &p1, Vector3f const &p2, int lineColor);
	void draw_triangle(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, int color, bool transparent);
	void draw_ground();
	void render_opaque		(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, bool textured, int tc, Vector2f const &tx1, Vector2f const &tx2, Vector2f const &tx3, Matrix22 const &txm);
	void render_transparent	(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, bool textured, int tc, Vector2f const &tx1, Vector2f const &tx2, Vector2f const &tx3, Matrix22 const &txm);
	void render();
	void show();
	void			rebuffer1			();
	void			rebuffer2			();
	void			rebuffer3			();
	void			rebuffer4			();
	void pushTexture(int *texture_);
	void enqueueTextures();
	int* popTexture();
	void clearTextures();
	void finish();
} nlsw;
void				NonlinearSW::initiate()
{
	texture=0;
	rgbn=h*w;
	wbuffer=(float*)malloc(rgbn*sizeof(float)), xbuffer=(float*)malloc(rgbn*sizeof(float)), ybuffer=(float*)malloc(rgbn*sizeof(float)), zbuffer=(float*)malloc(rgbn*sizeof(float));
	ghMemDC=CreateCompatibleDC(ghDC);
	tagBITMAPINFO bmpInfo={{sizeof(tagBITMAPINFOHEADER), w, -h, 1, 32, BI_RGB, 0, 0, 0, 0, 0}};
	hBitmap=(HBITMAP__*)SelectObject(ghMemDC, CreateDIBSection(0, &bmpInfo, DIB_RGB_COLORS, (void**)&rgb, 0, 0));
	SetBkMode(ghMemDC, TRANSPARENT);
}
void				NonlinearSW::resize()
{
	rgbn=h*w;
	wbuffer=(float*)realloc(wbuffer, rgbn*sizeof(float)), xbuffer=(float*)realloc(xbuffer, rgbn*sizeof(float)), ybuffer=(float*)realloc(ybuffer, rgbn*sizeof(float)), zbuffer=(float*)realloc(zbuffer, rgbn*sizeof(float));
	(this->*rebuffer)();
	DeleteObject((HBITMAP__*)SelectObject(ghMemDC, hBitmap));
	tagBITMAPINFO bmpInfo={{sizeof(tagBITMAPINFOHEADER), w, -h, 1, 32, BI_RGB, 0, 0, 0, 0, 0}};
	hBitmap=(HBITMAP__*)SelectObject(ghMemDC, CreateDIBSection(0, &bmpInfo, DIB_RGB_COLORS, (void**)&rgb, 0, 0));
	SetBkMode(ghMemDC, TRANSPARENT);
}
void				NonlinearSW::changeFov(){(this->*rebuffer)();}
void				NonlinearSW::clear_screen(){memset(rgb, 0xFF, rgbn*sizeof(int)), memset(wbuffer, 0, rgbn*sizeof(float));}
void				NonlinearSW::print(int x, int y, const char *a, ...){GUIVPrint(ghMemDC, x, y, a, (char*)(&a+1));}
void				NonlinearSW::print(Vector3f const &p, const char *a, ...)
{
	Vector2f s;
	Vector3f c;
	world_to_screen(p, s, c);
	if(c.z>0)
		GUIVPrint(ghMemDC, (int)s.x, (int)s.y, a, (char*)(&a+1));
}
void				NonlinearSW::draw_point(Vector3f const &p, int color){}
void				NonlinearSW::draw_line(float x1, float y1, float x2, float y2, int lineColor){}
void				NonlinearSW::draw_line(Vector3f const &p1, Vector3f const &p2, int lineColor)
{
	//Vector3f cp1, cp2;
	//world_to_camera(p1, cp1), world_to_camera(p2, cp2);
	//for(int ky=0;ky<h;++ky)
	//{
	//	for(int kx=0;kx<w;++kx)
	//	{
	//	}
	//}
}
void				NonlinearSW::draw_triangle(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, int color, bool transparent)
{
	Vector2f zero;
	if(transparent)
		render_transparent(p1, p2, p3, false, color, zero, zero, zero, Matrix22());
	else
		render_opaque(p1, p2, p3, false, color, zero, zero, zero, Matrix22());
}
void				NonlinearSW::draw_ground()
{
}
void				NonlinearSW::render_opaque		(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, bool textured, int tc, Vector2f const &tx1, Vector2f const &tx2, Vector2f const &tx3, Matrix22 const &txm)
{
	float cpt, v1[5], v2[5], v3[5], admittance;
//	float dx, dy, dz, cpt, v1[5], v2[5], v3[5], admittance;
	world_to_camera(p1, v1+2), world_to_camera(p2, v2+2), world_to_camera(p3, v3+2);
	//dx=p.p1.x-cam.x, dy=p.p1.y-cam.y, dz=p.p1.z-cam.z, cpt=dx*cax+dy*sax, v1[2]=dx*sax-dy*cax, v1[3]=cpt*say-dz*cay, v1[4]=cpt*cay+dz*say;
	//dx=p.p2.x-cam.x, dy=p.p2.y-cam.y, dz=p.p2.z-cam.z, cpt=dx*cax+dy*sax, v2[2]=dx*sax-dy*cax, v2[3]=cpt*say-dz*cay, v2[4]=cpt*cay+dz*say;
	//dx=p.p3.x-cam.x, dy=p.p3.y-cam.y, dz=p.p3.z-cam.z, cpt=dx*cax+dy*sax, v3[2]=dx*sax-dy*cax, v3[3]=cpt*say-dz*cay, v3[4]=cpt*cay+dz*say;
	float	au12=v2[2]-v1[2],			bu12=v2[3]-v1[3],			cu12=v2[4]-v1[4];									//u12	=<12>
	float	aux3=v3[2]-v1[2],			bux3=v3[3]-v1[3],			cux3=v3[4]-v1[4];									//ux3	=<13>
	float	a=bu12*cux3-bux3*cu12,		b=aux3*cu12-au12*cux3,		c=au12*bux3-aux3*bu12;								//abc	=<n>	=<12>x<13>
	float	t=a*v1[2]+b*v1[3]+c*v1[4], t2=t*t;
	if(!t)return;
	float	m12=sqrt(au12*au12+bu12*bu12+cu12*cu12);			au12/=m12,			bu12/=m12,			cu12/=m12;		//u12	=<u12>	=<12>/|12|
			cpt=aux3*au12+bux3*bu12+cux3*cu12;					aux3-=cpt*au12,		bux3-=cpt*bu12,		cux3-=cpt*cu12;	//ux3	=<x3>	=<13>-<u12>.<13><u12>
			cpt=sqrt(aux3*aux3+bux3*bux3+cux3*cux3);			aux3/=cpt,			bux3/=cpt,			cux3/=cpt;		//ux3	=<ux3>	=<x3>/|x3|
	float	ux3_1=aux3*v1[2]+bux3*v1[3]+cux3*v1[4];
	float	u12_1=au12*v1[2]+bu12*v1[3]+cu12*v1[4];
	float	Xp3=aux3*v3[2]+bux3*v3[3]+cux3*v3[4]-ux3_1;
	float	Yp3=au12*v3[2]+bu12*v3[3]+cu12*v3[4]-u12_1;
	float	A1=t*aux3-a*ux3_1, A2=t*bux3-b*ux3_1, A3=t*cux3-c*ux3_1;
	float	A4=t*au12-a*u12_1, A5=t*bu12-b*u12_1, A6=t*cu12-c*u12_1;
	if(textured)
	{
		float B1=txm.a*A1+txm.b*A4+tx1.x*a, B2=txm.a*A2+txm.b*A5+tx1.x*b, B3=txm.a*A3+txm.b*A6+tx1.x*c;
		float B4=txm.c*A1+txm.d*A4+tx1.y*a, B5=txm.c*A2+txm.d*A5+tx1.y*b, B6=txm.c*A3+txm.d*A6+tx1.y*c;
		float a12=tx1.y-tx2.y, b12=tx2.x-tx1.x, c12=-a12*tx1.x-b12*tx1.y, d12_3=a12*tx3.x+b12*tx3.y+c12;	a12/=d12_3, b12/=d12_3, c12/=d12_3;
		float a23=tx2.y-tx3.y, b23=tx3.x-tx2.x, c23=-a23*tx2.x-b23*tx2.y, d23_1=a23*tx1.x+b23*tx1.y+c23;	a23/=d23_1, b23/=d23_1, c23/=d23_1;
		float a31=tx3.y-tx1.y, b31=tx1.x-tx3.x, c31=-a31*tx3.x-b31*tx3.y, d31_2=a31*tx2.x+b31*tx2.y+c31;	a31/=d31_2, b31/=d31_2, c31/=d31_2;
		int *txk=texture[tc], txhk=txh[tc], txwk=txw[tc];
	//	b/=a, c/=a, B1/=a, B2/=a, B3/=a, B4/=a, B5/=a, B6/=a;//no: cpt
		if(SSE4_1)
		{
			int k2=0;
			__m128 s_txwk=_mm_set1_ps((float)txwk), s_txhk=_mm_set1_ps((float)txhk);
			__m128i i_txwk=_mm_set1_epi32(txwk), i_txhk=_mm_set1_epi32(txhk);
			__m128 m_a=_mm_set1_ps(a), m_b=_mm_set1_ps(b), m_c=_mm_set1_ps(c);
			__m128 m_t=_mm_set1_ps(t);
			__m128 m_B1=_mm_set1_ps(B1), m_B2=_mm_set1_ps(B2), m_B3=_mm_set1_ps(B3);
			__m128 m_B4=_mm_set1_ps(B4), m_B5=_mm_set1_ps(B5), m_B6=_mm_set1_ps(B6);
			__m128 m_a12=_mm_set1_ps(a12), m_b12=_mm_set1_ps(b12), m_c12=_mm_set1_ps(-c12);
			__m128 m_a23=_mm_set1_ps(a23), m_b23=_mm_set1_ps(b23), m_c23=_mm_set1_ps(-c23);
			__m128 m_a31=_mm_set1_ps(a31), m_b31=_mm_set1_ps(b31), m_c31=_mm_set1_ps(-c31);
			for(;k2+4<rgbn;k2+=4)
			{
				__m128 xb=_mm_loadu_ps(xbuffer+k2), yb=_mm_loadu_ps(ybuffer+k2), zb=_mm_loadu_ps(zbuffer+k2);
				__m128 den=_mm_mul_ps(m_a, xb);
				den=_mm_add_ps(den, _mm_mul_ps(m_b, yb));
				den=_mm_add_ps(den, _mm_mul_ps(m_c, zb));
				__m128 Xtx=_mm_mul_ps(m_B1, xb);
				__m128 Ytx=_mm_mul_ps(m_B4, xb);//y
				Xtx=_mm_add_ps(Xtx, _mm_mul_ps(m_B2, yb));
				Ytx=_mm_add_ps(Ytx, _mm_mul_ps(m_B5, yb));//y
				Xtx=_mm_add_ps(Xtx, _mm_mul_ps(m_B3, zb));
				Ytx=_mm_add_ps(Ytx, _mm_mul_ps(m_B6, zb));//y
				Xtx=_mm_div_ps(Xtx, den);
				Ytx=_mm_div_ps(Ytx, den);//y

				__m128 t1=_mm_mul_ps(den, m_t);
				__m128i cmp_t1=_mm_castps_si128(_mm_cmpge_ps(t1, _mm_setzero_ps()));

				__m128 l1=_mm_mul_ps(m_a12, Xtx);
				l1=_mm_add_ps(l1, _mm_mul_ps(m_b12, Ytx));
				__m128i cmp_l1=_mm_castps_si128(_mm_cmpgt_ps(l1, m_c12));

				__m128 l2=_mm_mul_ps(m_a23, Xtx);
				l2=_mm_add_ps(l2, _mm_mul_ps(m_b23, Ytx));
				__m128i cmp_l2=_mm_castps_si128(_mm_cmpgt_ps(l2, m_c23));

				__m128 l3=_mm_mul_ps(m_a31, Xtx);
				l3=_mm_add_ps(l3, _mm_mul_ps(m_b31, Ytx));
				__m128i cmp_l3=_mm_castps_si128(_mm_cmpgt_ps(l3, m_c31));

				__m128i cmp=_mm_and_si128(cmp_t1, cmp_l1);
				cmp=_mm_and_si128(cmp, cmp_l2);
				cmp=_mm_and_si128(cmp, cmp_l3);

				__m128i cmp2=_mm_shuffle_epi32(cmp, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
				__m128i cmp_or=_mm_or_si128(cmp, cmp2);//r3|r1, r2|r0, r1|r3, r0|r2
				cmp2=_mm_shuffle_epi32(cmp_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
				cmp_or=_mm_or_si128(cmp_or, cmp2);//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
				if(cmp_or.m128i_i32[0])
				{
					__m128 wb=_mm_loadu_ps(wbuffer+k2);
					__m128 m_adm=_mm_div_ps(den, m_t);
					__m128i c_a=_mm_castps_si128(_mm_cmpgt_ps(m_adm, wb));

					__m128i c_a2=_mm_shuffle_epi32(c_a, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
					__m128i c_a_or=_mm_or_si128(c_a, c_a2);//r3|r1, r2|r0, r1|r3, r0|r2
					c_a2=_mm_shuffle_epi32(c_a_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
					c_a_or=_mm_or_si128(c_a_or, c_a2);//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
					if(c_a_or.m128i_i32[0])
					{
						__m128i condition=_mm_and_si128(cmp, c_a);
						__m128i condition2=_mm_xor_si128(condition, _mm_set1_epi32(-1));

						wb=_mm_castsi128_ps(_mm_and_si128(condition2, _mm_castps_si128(wb)));
						__m128 wb2=_mm_castsi128_ps(_mm_and_si128(condition, _mm_castps_si128(m_adm)));
						wb2=_mm_castsi128_ps(_mm_or_si128(_mm_castps_si128(wb2), _mm_castps_si128(wb)));
						_mm_storeu_ps(wbuffer+k2, wb2);

						//mod operation		x%q = x-q*floor(x/q)
						__m128 Xtx_w=_mm_div_ps(Xtx, s_txwk);
						__m128 Ytx_h=_mm_div_ps(Ytx, s_txhk);

						Xtx_w=_mm_floor_ps(Xtx_w);//floor(x/txw)	SSE4.1
						Ytx_h=_mm_floor_ps(Ytx_h);

						Xtx_w=_mm_mul_ps(Xtx_w, s_txwk);//floor(x/txw)*txw
						Ytx_h=_mm_mul_ps(Ytx_h, s_txhk);

						Xtx=_mm_sub_ps(Xtx, Xtx_w);//x-floor(x/txw)*txw
						Ytx=_mm_sub_ps(Ytx, Ytx_h);
						
						Xtx=_mm_floor_ps(Xtx);//99.9 -> 99 not 100	SSE4.1
						Ytx=_mm_floor_ps(Ytx);

						__m128i i_Xtx=_mm_cvtps_epi32(Xtx);
						__m128i i_Ytx=_mm_cvtps_epi32(Ytx);

						__m128i tx_idx=_mm_mullo_epi32(i_txwk, i_Ytx);//SSE4.1
						tx_idx=_mm_add_epi32(tx_idx, i_Xtx);

						__m128i m_tx=_mm_set_epi32(txk[tx_idx.m128i_i32[3]], txk[tx_idx.m128i_i32[2]], txk[tx_idx.m128i_i32[1]], txk[tx_idx.m128i_i32[0]]);
						__m128i m_rgb=_mm_loadu_si128((__m128i*)(rgb+k2));
						m_tx=_mm_and_si128(m_tx, condition);
						m_rgb=_mm_and_si128(m_rgb, condition2);
						m_tx=_mm_or_si128(m_tx, m_rgb);
						_mm_storeu_si128((__m128i*)(rgb+k2), m_tx);
					}
				}
			}
			for(;k2<rgbn;++k2)
			{
				cpt=a*xbuffer[k2]+b*ybuffer[k2]+c*zbuffer[k2];//398ms fullscreen
				float Xtx=(B1*xbuffer[k2]+B2*ybuffer[k2]+B3*zbuffer[k2])/cpt;
				float Ytx=(B4*xbuffer[k2]+B5*ybuffer[k2]+B6*zbuffer[k2])/cpt;
				if(cpt*t>=0&&a12*Xtx+b12*Ytx+c12>0&&a23*Xtx+b23*Ytx+c23>0&&a31*Xtx+b31*Ytx+c31>0)
				{
					admittance=cpt/t;
					if(admittance>wbuffer[k2])
					{
						int Xtx2=int(Xtx)%txwk;if(Xtx2<0)Xtx2+=txwk;
						int Ytx2=int(Ytx)%txhk;if(Ytx2<0)Ytx2+=txhk;
						wbuffer[k2]=admittance, rgb[k2]=txk[Xtx2+Ytx2*txwk];
					}
				}
			}
		}
		else//textured, no SSE4.1
		{
			const __m128 half=_mm_set1_ps(0.5f);
			const __m128i mask2_lo=_mm_set_epi32(0, -1, 0, -1), mask2_hi=_mm_set_epi32(-1, 0, -1, 0);
			__m128 s_txwk=_mm_set1_ps((float)txwk), s_txhk=_mm_set1_ps((float)txhk);
			__m128i i_txwk=_mm_set1_epi32(txwk), i_txhk=_mm_set1_epi32(txhk);
			__m128 m_a=_mm_set1_ps(a), m_b=_mm_set1_ps(b), m_c=_mm_set1_ps(c);
			__m128 m_t=_mm_set1_ps(t);
			__m128 m_B1=_mm_set1_ps(B1), m_B2=_mm_set1_ps(B2), m_B3=_mm_set1_ps(B3);
			__m128 m_B4=_mm_set1_ps(B4), m_B5=_mm_set1_ps(B5), m_B6=_mm_set1_ps(B6);
			__m128 m_a12=_mm_set1_ps(a12), m_b12=_mm_set1_ps(b12), m_c12=_mm_set1_ps(-c12);
			__m128 m_a23=_mm_set1_ps(a23), m_b23=_mm_set1_ps(b23), m_c23=_mm_set1_ps(-c23);
			__m128 m_a31=_mm_set1_ps(a31), m_b31=_mm_set1_ps(b31), m_c31=_mm_set1_ps(-c31);
			int k2=0;
			for(;k2+4<rgbn;k2+=4)
			{
				__m128 xb=_mm_loadu_ps(xbuffer+k2), yb=_mm_loadu_ps(ybuffer+k2), zb=_mm_loadu_ps(zbuffer+k2);
				__m128 den=_mm_mul_ps(m_a, xb);
				den=_mm_add_ps(den, _mm_mul_ps(m_b, yb));
				den=_mm_add_ps(den, _mm_mul_ps(m_c, zb));
				__m128 Xtx=_mm_mul_ps(m_B1, xb);
				__m128 Ytx=_mm_mul_ps(m_B4, xb);//y
				Xtx=_mm_add_ps(Xtx, _mm_mul_ps(m_B2, yb));
				Ytx=_mm_add_ps(Ytx, _mm_mul_ps(m_B5, yb));//y
				Xtx=_mm_add_ps(Xtx, _mm_mul_ps(m_B3, zb));
				Ytx=_mm_add_ps(Ytx, _mm_mul_ps(m_B6, zb));//y
				Xtx=_mm_div_ps(Xtx, den);
				Ytx=_mm_div_ps(Ytx, den);//y

				__m128 t1=_mm_mul_ps(den, m_t);
				__m128i cmp_t1=_mm_castps_si128(_mm_cmpge_ps(t1, _mm_setzero_ps()));

				__m128 l1=_mm_mul_ps(m_a12, Xtx);
				l1=_mm_add_ps(l1, _mm_mul_ps(m_b12, Ytx));
				__m128i cmp_l1=_mm_castps_si128(_mm_cmpgt_ps(l1, m_c12));

				__m128 l2=_mm_mul_ps(m_a23, Xtx);
				l2=_mm_add_ps(l2, _mm_mul_ps(m_b23, Ytx));
				__m128i cmp_l2=_mm_castps_si128(_mm_cmpgt_ps(l2, m_c23));

				__m128 l3=_mm_mul_ps(m_a31, Xtx);
				l3=_mm_add_ps(l3, _mm_mul_ps(m_b31, Ytx));
				__m128i cmp_l3=_mm_castps_si128(_mm_cmpgt_ps(l3, m_c31));

				__m128i cmp=_mm_and_si128(cmp_t1, cmp_l1);
				cmp=_mm_and_si128(cmp, cmp_l2);
				cmp=_mm_and_si128(cmp, cmp_l3);

				__m128i cmp2=_mm_shuffle_epi32(cmp, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
				__m128i cmp_or=_mm_or_si128(cmp, cmp2);//r3|r1, r2|r0, r1|r3, r0|r2
				cmp2=_mm_shuffle_epi32(cmp_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
				cmp_or=_mm_or_si128(cmp_or, cmp2);//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
				if(cmp_or.m128i_i32[0])
				{
					__m128 wb=_mm_loadu_ps(wbuffer+k2);
					__m128 m_adm=_mm_div_ps(den, m_t);
					__m128i c_a=_mm_castps_si128(_mm_cmpgt_ps(m_adm, wb));

					__m128i c_a2=_mm_shuffle_epi32(c_a, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
					__m128i c_a_or=_mm_or_si128(c_a, c_a2);//r3|r1, r2|r0, r1|r3, r0|r2
					c_a2=_mm_shuffle_epi32(c_a_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
					c_a_or=_mm_or_si128(c_a_or, c_a2);//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
					if(c_a_or.m128i_i32[0])
					{
						__m128i condition=_mm_and_si128(cmp, c_a);
						__m128i condition2=_mm_xor_si128(condition, _mm_set1_epi32(-1));

						wb=_mm_castsi128_ps(_mm_and_si128(condition2, _mm_castps_si128(wb)));
						__m128 wb2=_mm_castsi128_ps(_mm_and_si128(condition, _mm_castps_si128(m_adm)));
						wb2=_mm_castsi128_ps(_mm_or_si128(_mm_castps_si128(wb2), _mm_castps_si128(wb)));
						_mm_storeu_ps(wbuffer+k2, wb2);

						//mod operation		x%q = x-q*floor(x/q)
						__m128 Xtx_w=_mm_div_ps(Xtx, s_txwk);
						__m128 Ytx_h=_mm_div_ps(Ytx, s_txhk);

						//Xtx_w=_mm_floor_ps(Xtx_w);//floor(x/txw)	SSE4.1
						//Ytx_h=_mm_floor_ps(Ytx_h);
						Xtx_w=_mm_sub_ps(Xtx_w, half);//floor(x/q)
						Ytx_h=_mm_sub_ps(Ytx_h, half);
						__m128i i_Xtx_w=_mm_cvtps_epi32(Xtx_w);
						__m128i i_Ytx_h=_mm_cvtps_epi32(Ytx_h);
						Xtx_w=_mm_cvtepi32_ps(i_Xtx_w);
						Ytx_h=_mm_cvtepi32_ps(i_Ytx_h);

						Xtx_w=_mm_mul_ps(Xtx_w, s_txwk);//floor(x/txw)*txw
						Ytx_h=_mm_mul_ps(Ytx_h, s_txhk);

						Xtx=_mm_sub_ps(Xtx, Xtx_w);//x-floor(x/txw)*txw
						Ytx=_mm_sub_ps(Ytx, Ytx_h);
						
						//Xtx=_mm_floor_ps(Xtx);//99.9 -> 99 not 100	SSE4.1
						//Ytx=_mm_floor_ps(Ytx);
						Xtx=_mm_sub_ps(Xtx, half);
						Ytx=_mm_sub_ps(Ytx, half);

						__m128i i_Xtx=_mm_cvtps_epi32(Xtx);
						__m128i i_Ytx=_mm_cvtps_epi32(Ytx);

						//__m128i tx_idx=_mm_mullo_epi32(i_txwk, i_Ytx);//SSE4.1
						__m128i tx_idx_lo=_mm_mul_epu32(i_Ytx, i_txwk);
						i_Ytx=_mm_shuffle_epi32(i_Ytx, _MM_SHUFFLE(2, 3, 0, 1));//y2 y3 y0 y1
						__m128i tx_idx_hi=_mm_mul_epu32(i_Ytx, i_txwk);
						tx_idx_hi=_mm_shuffle_epi32(tx_idx_hi, _MM_SHUFFLE(2, 3, 0, 1));
						tx_idx_lo=_mm_and_si128(tx_idx_lo, mask2_lo);
						tx_idx_hi=_mm_and_si128(tx_idx_hi, mask2_hi);
						__m128i tx_idx=_mm_or_si128(tx_idx_lo, tx_idx_hi);

						tx_idx=_mm_add_epi32(tx_idx, i_Xtx);

						__m128i m_tx=_mm_set_epi32(txk[tx_idx.m128i_i32[3]], txk[tx_idx.m128i_i32[2]], txk[tx_idx.m128i_i32[1]], txk[tx_idx.m128i_i32[0]]);
						__m128i m_rgb=_mm_loadu_si128((__m128i*)(rgb+k2));
						m_tx=_mm_and_si128(m_tx, condition);
						m_rgb=_mm_and_si128(m_rgb, condition2);
						m_tx=_mm_or_si128(m_tx, m_rgb);
						_mm_storeu_si128((__m128i*)(rgb+k2), m_tx);
					}
				}
			}
			for(;k2<rgbn;++k2)
			{
				cpt=a*xbuffer[k2]+b*ybuffer[k2]+c*zbuffer[k2];//398ms fullscreen
				float Xtx=(B1*xbuffer[k2]+B2*ybuffer[k2]+B3*zbuffer[k2])/cpt;
				float Ytx=(B4*xbuffer[k2]+B5*ybuffer[k2]+B6*zbuffer[k2])/cpt;
				if(cpt*t>=0&&a12*Xtx+b12*Ytx+c12>0&&a23*Xtx+b23*Ytx+c23>0&&a31*Xtx+b31*Ytx+c31>0)
				{
					admittance=cpt/t;
					if(admittance>wbuffer[k2])
					{
						int Xtx2=int(Xtx)%txwk;if(Xtx2<0)Xtx2+=txwk;
						int Ytx2=int(Ytx)%txhk;if(Ytx2<0)Ytx2+=txhk;
						wbuffer[k2]=admittance, rgb[k2]=txk[Xtx2+Ytx2*txwk];
					}
				}
			}
		/*	for(int k2=0;k2<rgbn;++k2)
			{
			//	cpt=xbuffer[k2]+b*ybuffer[k2]+c*zbuffer[k2];//no: cpt
				//auto &ux=xbuffer[k2], &uy=ybuffer[k2], &uz=zbuffer[k2];//402ms fullscreen
				//cpt=a*ux+b*uy+c*zbuffer[k2];
				//float Xtx=(B1*ux+B2*uy+B3*uz)/cpt;
				//float Ytx=(B4*ux+B5*uy+B6*uz)/cpt;
				cpt=a*xbuffer[k2]+b*ybuffer[k2]+c*zbuffer[k2];//398ms fullscreen
				float Xtx=(B1*xbuffer[k2]+B2*ybuffer[k2]+B3*zbuffer[k2])/cpt;
				float Ytx=(B4*xbuffer[k2]+B5*ybuffer[k2]+B6*zbuffer[k2])/cpt;
				if(cpt*t>=0&&a12*Xtx+b12*Ytx+c12>0&&a23*Xtx+b23*Ytx+c23>0&&a31*Xtx+b31*Ytx+c31>0)
				{
					admittance=cpt/t;
					if(admittance>wbuffer[k2])
					{
						int Xtx2=int(Xtx)%txwk;if(Xtx2<0)Xtx2+=txwk;
						int Ytx2=int(Ytx)%txhk;if(Ytx2<0)Ytx2+=txhk;
						wbuffer[k2]=admittance, rgb[k2]=txk[Xtx2+Ytx2*txwk];
					}
				}
				//admittance=cpt/t;
				//if(admittance>=0&&a12*Xtx+b12*Ytx+c12>0&&a23*Xtx+b23*Ytx+c23>0&&a31*Xtx+b31*Ytx+c31>0&&admittance>wbuffer[k2])
				//{
				//	int Xtx2=int(Xtx)%txwk;if(Xtx2<0)Xtx2+=txwk;
				//	int Ytx2=int(Ytx)%txhk;if(Ytx2<0)Ytx2+=txhk;
				//	wbuffer[k2]=admittance, rgb[k2]=txk[Xtx2+Ytx2*txwk];
				//}
				//if(t*cpt>=0&&a12*Xtx+b12*Ytx+c12>0&&a23*Xtx+b23*Ytx+c23>0&&a31*Xtx+b31*Ytx+c31>0)
				//{
				//	admittance=cpt*cpt/t2;
				//	if(admittance>wbuffer[k2])
				//	{
				//		int Xtx2=int(Xtx)%txwk;if(Xtx2<0)Xtx2+=txwk;
				//		int Ytx2=int(Ytx)%txhk;if(Ytx2<0)Ytx2+=txhk;
				//		wbuffer[k2]=admittance, rgb[k2]=txk[Xtx2+Ytx2*txwk];
				//	}
				//}
			}//*/
		}
	}
	else//solid color
	{
		float a12=m12,											r12=a12*Xp3;		a12/=r12;
		float a23=Yp3-m12,	b23=-Xp3,	c23=-b23*m12,			r23=c23;			a23/=r23, b23/=r23, c23/=r23;
		float a31=-Yp3,		b31=Xp3,	c31=-a31*Xp3-b31*Yp3,	r31=b31*m12+c31;	a31/=r31, b31/=r31, c31/=r31;
		int color=tc;
	//	b/=a, c/=a, A1/=a, A2/=a, A3/=a, A4/=a, A5/=a, A6/=a;//cpt
		//if(SSE4_1)
		//{
			int k2=0;
			__m128i m_color=_mm_set1_epi32(color);
			__m128 m_a=_mm_set1_ps(a), m_b=_mm_set1_ps(b), m_c=_mm_set1_ps(c);
			__m128 m_t=_mm_set1_ps(t);
			__m128 m_A1=_mm_set1_ps(A1), m_A2=_mm_set1_ps(A2), m_A3=_mm_set1_ps(A3);
			__m128 m_A4=_mm_set1_ps(A4), m_A5=_mm_set1_ps(A5), m_A6=_mm_set1_ps(A6);
			__m128 m_a12=_mm_set1_ps(a12);
			__m128 m_a23=_mm_set1_ps(a23), m_b23=_mm_set1_ps(b23), m_c23=_mm_set1_ps(-c23);
			__m128 m_a31=_mm_set1_ps(a31), m_b31=_mm_set1_ps(b31), m_c31=_mm_set1_ps(-c31);
			for(;k2+4<rgbn;k2+=4)
			{
				__m128 xb=_mm_loadu_ps(xbuffer+k2), yb=_mm_loadu_ps(ybuffer+k2), zb=_mm_loadu_ps(zbuffer+k2);
				__m128 den=_mm_mul_ps(m_a, xb);
				den=_mm_add_ps(den, _mm_mul_ps(m_b, yb));
				den=_mm_add_ps(den, _mm_mul_ps(m_c, zb));
				__m128 Xtx=_mm_mul_ps(m_A1, xb);
				__m128 Ytx=_mm_mul_ps(m_A4, xb);//y
				Xtx=_mm_add_ps(Xtx, _mm_mul_ps(m_A2, yb));
				Ytx=_mm_add_ps(Ytx, _mm_mul_ps(m_A5, yb));//y
				Xtx=_mm_add_ps(Xtx, _mm_mul_ps(m_A3, zb));
				Ytx=_mm_add_ps(Ytx, _mm_mul_ps(m_A6, zb));//y
				Xtx=_mm_div_ps(Xtx, den);
				Ytx=_mm_div_ps(Ytx, den);//y

				__m128 t1=_mm_mul_ps(den, m_t);
				__m128i cmp_t1=_mm_castps_si128(_mm_cmpge_ps(t1, _mm_setzero_ps()));

				__m128 l1=_mm_mul_ps(m_a12, Xtx);
				__m128i cmp_l1=_mm_castps_si128(_mm_cmpgt_ps(l1, _mm_setzero_ps()));

				__m128 l2=_mm_mul_ps(m_a23, Xtx);
				l2=_mm_add_ps(l2, _mm_mul_ps(m_b23, Ytx));
				__m128i cmp_l2=_mm_castps_si128(_mm_cmpgt_ps(l2, m_c23));

				__m128 l3=_mm_mul_ps(m_a31, Xtx);
				l3=_mm_add_ps(l3, _mm_mul_ps(m_b31, Ytx));
				__m128i cmp_l3=_mm_castps_si128(_mm_cmpgt_ps(l3, m_c31));

				__m128i cmp=_mm_and_si128(cmp_t1, cmp_l1);
				cmp=_mm_and_si128(cmp, cmp_l2);
				cmp=_mm_and_si128(cmp, cmp_l3);

				__m128i cmp2=_mm_shuffle_epi32(cmp, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
				__m128i cmp_or=_mm_or_si128(cmp, cmp2);//r3|r1, r2|r0, r1|r3, r0|r2
				cmp2=_mm_shuffle_epi32(cmp_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
				cmp_or=_mm_or_si128(cmp_or, cmp2);//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
				if(cmp_or.m128i_i32[0])
				{
					__m128 wb=_mm_loadu_ps(wbuffer+k2);
					__m128 m_adm=_mm_div_ps(den, m_t);
					__m128i c_a=_mm_castps_si128(_mm_cmpgt_ps(m_adm, wb));
					
					__m128i c_a2=_mm_shuffle_epi32(c_a, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
					__m128i c_a_or=_mm_or_si128(c_a, c_a2);//r3|r1, r2|r0, r1|r3, r0|r2
					c_a2=_mm_shuffle_epi32(c_a_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
					c_a_or=_mm_or_si128(c_a_or, c_a2);//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
					if(c_a_or.m128i_i32[0])
					{
						__m128i condition=_mm_and_si128(cmp, c_a);
						__m128i condition2=_mm_xor_si128(condition, _mm_set1_epi32(-1));

						wb=_mm_castsi128_ps(_mm_and_si128(condition2, _mm_castps_si128(wb)));
						__m128 wb2=_mm_castsi128_ps(_mm_and_si128(condition, _mm_castps_si128(m_adm)));
						wb2=_mm_castsi128_ps(_mm_or_si128(_mm_castps_si128(wb2), _mm_castps_si128(wb)));
						_mm_storeu_ps(wbuffer+k2, wb2);
						
						__m128i m_rgb=_mm_loadu_si128((__m128i*)(rgb+k2));
						__m128i m_color2=_mm_and_si128(m_color, condition);
						m_rgb=_mm_and_si128(m_rgb, condition2);
						m_color2=_mm_or_si128(m_color2, m_rgb);
						_mm_storeu_si128((__m128i*)(rgb+k2), m_color2);
					}
				}
			}
			for(;k2<rgbn;++k2)
			{
				cpt=a*xbuffer[k2]+b*ybuffer[k2]+c*zbuffer[k2];
				float Xp=(A1*xbuffer[k2]+A2*ybuffer[k2]+A3*zbuffer[k2])/cpt;
				float Yp=(A4*xbuffer[k2]+A5*ybuffer[k2]+A6*zbuffer[k2])/cpt;
				if(cpt*t>=0&&a12*Xp>0&&a23*Xp+b23*Yp+c23>0&&a31*Xp+b31*Yp+c31>0)
				{
					admittance=cpt/t;
					if(admittance>wbuffer[k2])
						wbuffer[k2]=admittance, rgb[k2]=color;
				}
			}
		//}
		//else
		//{
		//	for(int k2=0;k2<rgbn;++k2)
		//	{
		//	//	cpt=xbuffer[k2]+b*ybuffer[k2]+c*zbuffer[k2];//
		//		//auto &ux=xbuffer[k2], &uy=ybuffer[k2], &uz=zbuffer[k2];
		//		//cpt=a*ux+b*uy+c*zbuffer[k2];
		//		//float Xp=(A1*ux+A2*uy+A3*uz)/cpt;
		//		//float Yp=(A4*ux+A5*uy+A6*uz)/cpt;
		//		cpt=a*xbuffer[k2]+b*ybuffer[k2]+c*zbuffer[k2];
		//		float Xp=(A1*xbuffer[k2]+A2*ybuffer[k2]+A3*zbuffer[k2])/cpt;
		//		float Yp=(A4*xbuffer[k2]+A5*ybuffer[k2]+A6*zbuffer[k2])/cpt;
		//		if(cpt*t>=0&&a12*Xp>0&&a23*Xp+b23*Yp+c23>0&&a31*Xp+b31*Yp+c31>0)
		//		{
		//			admittance=cpt/t;
		//			if(admittance>wbuffer[k2])
		//				wbuffer[k2]=admittance, rgb[k2]=color;
		//		}
		//		//admittance=cpt/t;
		//		//if(admittance>=0&&a12*Xp>0&&a23*Xp+b23*Yp+c23>0&&a31*Xp+b31*Yp+c31>0&&admittance>wbuffer[k2])
		//		//	wbuffer[k2]=admittance, rgb[k2]=color;
		//		//if(t*cpt>=0&&a12*Xp>0&&a23*Xp+b23*Yp+c23>0&&a31*Xp+b31*Yp+c31>0)
		//		//{
		//		//	admittance=cpt*cpt/t2;
		//		//	if(admittance>wbuffer[k2])
		//		//		wbuffer[k2]=admittance, rgb[k2]=color;
		//		//}
		//	}
		//}
	}
}
void				NonlinearSW::render_transparent	(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, bool textured, int tc, Vector2f const &tx1, Vector2f const &tx2, Vector2f const &tx3, Matrix22 const &txm)
{
	float cpt, v1[5], v2[5], v3[5], admittance;
//	float dx, dy, dz, cpt, v1[5], v2[5], v3[5], admittance;
	world_to_camera(p1, v1+2), world_to_camera(p2, v2+2), world_to_camera(p3, v3+2);
	float	au12=v2[2]-v1[2],			bu12=v2[3]-v1[3],			cu12=v2[4]-v1[4];									//u12	=<12>
	float	aux3=v3[2]-v1[2],			bux3=v3[3]-v1[3],			cux3=v3[4]-v1[4];									//ux3	=<13>
	float	a=bu12*cux3-bux3*cu12,		b=aux3*cu12-au12*cux3,		c=au12*bux3-aux3*bu12;								//abc	=<n>	=<12>x<13>
	float	t=a*v1[2]+b*v1[3]+c*v1[4], t2=t*t;
	if(!t)return;
	float	m12=sqrt(au12*au12+bu12*bu12+cu12*cu12);			au12/=m12,			bu12/=m12,			cu12/=m12;		//u12	=<u12>	=<12>/|12|
			cpt=aux3*au12+bux3*bu12+cux3*cu12;					aux3-=cpt*au12,		bux3-=cpt*bu12,		cux3-=cpt*cu12;	//ux3	=<x3>	=<13>-<u12>.<13><u12>
			cpt=sqrt(aux3*aux3+bux3*bux3+cux3*cux3);			aux3/=cpt,			bux3/=cpt,			cux3/=cpt;		//ux3	=<ux3>	=<x3>/|x3|
	float	ux3_1=aux3*v1[2]+bux3*v1[3]+cux3*v1[4];
	float	u12_1=au12*v1[2]+bu12*v1[3]+cu12*v1[4];
	float	Xp3=aux3*v3[2]+bux3*v3[3]+cux3*v3[4]-ux3_1;
	float	Yp3=au12*v3[2]+bu12*v3[3]+cu12*v3[4]-u12_1;
	float	A1=t*aux3-a*ux3_1, A2=t*bux3-b*ux3_1, A3=t*cux3-c*ux3_1;
	float	A4=t*au12-a*u12_1, A5=t*bu12-b*u12_1, A6=t*cu12-c*u12_1;
	if(textured)
	{
		float B1=txm.a*A1+txm.b*A4+tx1.x*a, B2=txm.a*A2+txm.b*A5+tx1.x*b, B3=txm.a*A3+txm.b*A6+tx1.x*c;
		float B4=txm.c*A1+txm.d*A4+tx1.y*a, B5=txm.c*A2+txm.d*A5+tx1.y*b, B6=txm.c*A3+txm.d*A6+tx1.y*c;
		float a12=tx1.y-tx2.y, b12=tx2.x-tx1.x, c12=-a12*tx1.x-b12*tx1.y, d12_3=a12*tx3.x+b12*tx3.y+c12;	a12/=d12_3, b12/=d12_3, c12/=d12_3;
		float a23=tx2.y-tx3.y, b23=tx3.x-tx2.x, c23=-a23*tx2.x-b23*tx2.y, d23_1=a23*tx1.x+b23*tx1.y+c23;	a23/=d23_1, b23/=d23_1, c23/=d23_1;
		float a31=tx3.y-tx1.y, b31=tx1.x-tx3.x, c31=-a31*tx3.x-b31*tx3.y, d31_2=a31*tx2.x+b31*tx2.y+c31;	a31/=d31_2, b31/=d31_2, c31/=d31_2;
		int *txk=texture[tc], txhk=txh[tc], txwk=txw[tc];
	//	b/=a, c/=a, B1/=a, B2/=a, B3/=a, B4/=a, B5/=a, B6/=a;//cpt
		if(SSE4_1)
		{
			int k2=0;
			__m128i const shuffle_lo=_mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 15, 13, 11, 9, 7, 5, 3, 1);
			__m128i const shuffle_hi=_mm_set_epi8(15, 13, 11, 9, 7, 5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0);
			__m128i const mask_lo=_mm_set_epi32(0, 0, -1, -1);
			__m128i const mask_hi=_mm_set_epi32(-1, -1, 0, 0);
			__m128 s_txwk=_mm_set1_ps((float)txwk), s_txhk=_mm_set1_ps((float)txhk);
			__m128i i_txwk=_mm_set1_epi32(txwk), i_txhk=_mm_set1_epi32(txhk);
			__m128 m_a=_mm_set1_ps(a), m_b=_mm_set1_ps(b), m_c=_mm_set1_ps(c);
			__m128 m_t=_mm_set1_ps(t);
			__m128 m_B1=_mm_set1_ps(B1), m_B2=_mm_set1_ps(B2), m_B3=_mm_set1_ps(B3);
			__m128 m_B4=_mm_set1_ps(B4), m_B5=_mm_set1_ps(B5), m_B6=_mm_set1_ps(B6);
			__m128 m_a12=_mm_set1_ps(a12), m_b12=_mm_set1_ps(b12), m_c12=_mm_set1_ps(-c12);
			__m128 m_a23=_mm_set1_ps(a23), m_b23=_mm_set1_ps(b23), m_c23=_mm_set1_ps(-c23);
			__m128 m_a31=_mm_set1_ps(a31), m_b31=_mm_set1_ps(b31), m_c31=_mm_set1_ps(-c31);
			for(;k2+4<rgbn;k2+=4)
			{
				__m128 xb=_mm_loadu_ps(xbuffer+k2), yb=_mm_loadu_ps(ybuffer+k2), zb=_mm_loadu_ps(zbuffer+k2);
				__m128 den=_mm_mul_ps(m_a, xb);
				den=_mm_add_ps(den, _mm_mul_ps(m_b, yb));
				den=_mm_add_ps(den, _mm_mul_ps(m_c, zb));
				__m128 Xtx=_mm_mul_ps(m_B1, xb);
				__m128 Ytx=_mm_mul_ps(m_B4, xb);//y
				Xtx=_mm_add_ps(Xtx, _mm_mul_ps(m_B2, yb));
				Ytx=_mm_add_ps(Ytx, _mm_mul_ps(m_B5, yb));//y
				Xtx=_mm_add_ps(Xtx, _mm_mul_ps(m_B3, zb));
				Ytx=_mm_add_ps(Ytx, _mm_mul_ps(m_B6, zb));//y
				Xtx=_mm_div_ps(Xtx, den);
				Ytx=_mm_div_ps(Ytx, den);//y

				__m128 t1=_mm_mul_ps(den, m_t);
				__m128i cmp_t1=_mm_castps_si128(_mm_cmpge_ps(t1, _mm_setzero_ps()));

				__m128 l1=_mm_mul_ps(m_a12, Xtx);
				l1=_mm_add_ps(l1, _mm_mul_ps(m_b12, Ytx));
				__m128i cmp_l1=_mm_castps_si128(_mm_cmpgt_ps(l1, m_c12));

				__m128 l2=_mm_mul_ps(m_a23, Xtx);
				l2=_mm_add_ps(l2, _mm_mul_ps(m_b23, Ytx));
				__m128i cmp_l2=_mm_castps_si128(_mm_cmpgt_ps(l2, m_c23));

				__m128 l3=_mm_mul_ps(m_a31, Xtx);
				l3=_mm_add_ps(l3, _mm_mul_ps(m_b31, Ytx));
				__m128i cmp_l3=_mm_castps_si128(_mm_cmpgt_ps(l3, m_c31));

				__m128i cmp=_mm_and_si128(cmp_t1, cmp_l1);
				cmp=_mm_and_si128(cmp, cmp_l2);
				cmp=_mm_and_si128(cmp, cmp_l3);

				__m128i cmp2=_mm_shuffle_epi32(cmp, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
				__m128i cmp_or=_mm_or_si128(cmp, cmp2);//r3|r1, r2|r0, r1|r3, r0|r2
				cmp2=_mm_shuffle_epi32(cmp_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
				cmp_or=_mm_or_si128(cmp_or, cmp2);//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
				if(cmp_or.m128i_i32[0])
				{
					__m128 wb=_mm_loadu_ps(wbuffer+k2);
					__m128 m_adm=_mm_div_ps(den, m_t);
					__m128i c_a=_mm_castps_si128(_mm_cmpgt_ps(m_adm, wb));

					__m128i c_a2=_mm_shuffle_epi32(c_a, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
					__m128i c_a_or=_mm_or_si128(c_a, c_a2);//r3|r1, r2|r0, r1|r3, r0|r2
					c_a2=_mm_shuffle_epi32(c_a_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
					c_a_or=_mm_or_si128(c_a_or, c_a2);//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
					if(c_a_or.m128i_i32[0])
					{
						__m128i condition=_mm_and_si128(cmp, c_a);
						__m128i condition2=_mm_xor_si128(condition, _mm_set1_epi32(-1));

						wb=_mm_castsi128_ps(_mm_and_si128(condition2, _mm_castps_si128(wb)));
						__m128 wb2=_mm_castsi128_ps(_mm_and_si128(condition, _mm_castps_si128(m_adm)));
						wb2=_mm_castsi128_ps(_mm_or_si128(_mm_castps_si128(wb2), _mm_castps_si128(wb)));
						_mm_storeu_ps(wbuffer+k2, wb2);

						//mod operation		x%q = x-q*floor(x/q)
						__m128 Xtx_w=_mm_div_ps(Xtx, s_txwk);
						__m128 Ytx_h=_mm_div_ps(Ytx, s_txhk);

						Xtx_w=_mm_floor_ps(Xtx_w);//floor(x/txw)	SSE4.1
						Ytx_h=_mm_floor_ps(Ytx_h);

						Xtx_w=_mm_mul_ps(Xtx_w, s_txwk);//floor(x/txw)*txw
						Ytx_h=_mm_mul_ps(Ytx_h, s_txhk);

						Xtx=_mm_sub_ps(Xtx, Xtx_w);//x-floor(x/txw)*txw
						Ytx=_mm_sub_ps(Ytx, Ytx_h);
						
						Xtx=_mm_floor_ps(Xtx);//99.9 -> 99 not 100	SSE4.1
						Ytx=_mm_floor_ps(Ytx);

						__m128i i_Xtx=_mm_cvtps_epi32(Xtx);
						__m128i i_Ytx=_mm_cvtps_epi32(Ytx);

						__m128i tx_idx=_mm_mullo_epi32(i_txwk, i_Ytx);//SSE4.1
						tx_idx=_mm_add_epi32(tx_idx, i_Xtx);

						__m128i m_tx=_mm_set_epi32(txk[tx_idx.m128i_i32[3]], txk[tx_idx.m128i_i32[2]], txk[tx_idx.m128i_i32[1]], txk[tx_idx.m128i_i32[0]]);
						__m128i m_rgb=_mm_loadu_si128((__m128i*)(rgb+k2));
						__m128i m_rgb2=_mm_and_si128(m_rgb, condition2);
						m_tx=_mm_and_si128(m_tx, condition);
						m_rgb=_mm_and_si128(m_rgb, condition);

						__m128i const zero=_mm_setzero_si128();
						__m128i m_tx_lo=_mm_unpacklo_epi8(m_tx, zero);//{0, tx[7], ..., 0, tx[0]}
						__m128i m_tx_hi=_mm_unpackhi_epi8(m_tx, zero);//{0, tx[15], ..., 0, tx[8]}
						__m128i m_rgb_lo=_mm_unpacklo_epi8(m_rgb, zero);
						__m128i m_rgb_hi=_mm_unpackhi_epi8(m_rgb, zero);
						m_tx_lo=_mm_mullo_epi16(m_tx_lo, m_rgb_lo);//{{0, tx[7]}*{0, rgb[7]}, ..., {0, tx[0]}*{0, rgb[0]}}
						m_tx_hi=_mm_mullo_epi16(m_tx_hi, m_rgb_hi);
						m_tx_lo=_mm_shuffle_epi8(m_tx_lo, shuffle_lo);
						m_tx_hi=_mm_shuffle_epi8(m_tx_hi, shuffle_hi);
						m_tx_lo=_mm_and_si128(m_tx_lo, mask_lo);
						m_tx_hi=_mm_and_si128(m_tx_hi, mask_hi);
						m_tx=_mm_or_si128(m_tx_lo, m_tx_hi);

						m_tx=_mm_or_si128(m_tx, m_rgb2);
						_mm_storeu_si128((__m128i*)(rgb+k2), m_tx);
					}
				}
			}
			for(;k2<rgbn;++k2)
			{
				cpt=a*xbuffer[k2]+b*ybuffer[k2]+c*zbuffer[k2];
				float Xtx=(B1*xbuffer[k2]+B2*ybuffer[k2]+B3*zbuffer[k2])/cpt;
				float Ytx=(B4*xbuffer[k2]+B5*ybuffer[k2]+B6*zbuffer[k2])/cpt;
				if(cpt*t>=0&&a12*Xtx+b12*Ytx+c12>0&&a23*Xtx+b23*Ytx+c23>0&&a31*Xtx+b31*Ytx+c31>0)
				{
					admittance=cpt/t;
					if(admittance>wbuffer[k2])
					{
						int Xtx2=int(Xtx)%txwk;if(Xtx2<0)Xtx2+=txwk;
						int Ytx2=int(Ytx)%txhk;if(Ytx2<0)Ytx2+=txhk;
						((unsigned char*)&rgb[k2])[0]=((unsigned char*)&rgb[k2])[0]*((unsigned char*)&txk[Xtx2+Ytx2*txwk])[0]>>8;//b	//f(x,y)=x*y
						((unsigned char*)&rgb[k2])[1]=((unsigned char*)&rgb[k2])[1]*((unsigned char*)&txk[Xtx2+Ytx2*txwk])[1]>>8;//g
						((unsigned char*)&rgb[k2])[2]=((unsigned char*)&rgb[k2])[2]*((unsigned char*)&txk[Xtx2+Ytx2*txwk])[2]>>8;//r
					}
				}
			}
		}
		else//transparent textured, no SSE4.1
		{
			const __m128 half=_mm_set1_ps(0.5f);
			const __m128i mask2_lo=_mm_set_epi32(0, -1, 0, -1), mask2_hi=_mm_set_epi32(-1, 0, -1, 0);
			__m128i const shuffle_lo=_mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 15, 13, 11, 9, 7, 5, 3, 1);
			__m128i const shuffle_hi=_mm_set_epi8(15, 13, 11, 9, 7, 5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0);
			__m128i const mask_lo=_mm_set_epi32(0, 0, -1, -1);
			__m128i const mask_hi=_mm_set_epi32(-1, -1, 0, 0);
			__m128 s_txwk=_mm_set1_ps((float)txwk), s_txhk=_mm_set1_ps((float)txhk);
			__m128i i_txwk=_mm_set1_epi32(txwk), i_txhk=_mm_set1_epi32(txhk);
			__m128 m_a=_mm_set1_ps(a), m_b=_mm_set1_ps(b), m_c=_mm_set1_ps(c);
			__m128 m_t=_mm_set1_ps(t);
			__m128 m_B1=_mm_set1_ps(B1), m_B2=_mm_set1_ps(B2), m_B3=_mm_set1_ps(B3);
			__m128 m_B4=_mm_set1_ps(B4), m_B5=_mm_set1_ps(B5), m_B6=_mm_set1_ps(B6);
			__m128 m_a12=_mm_set1_ps(a12), m_b12=_mm_set1_ps(b12), m_c12=_mm_set1_ps(-c12);
			__m128 m_a23=_mm_set1_ps(a23), m_b23=_mm_set1_ps(b23), m_c23=_mm_set1_ps(-c23);
			__m128 m_a31=_mm_set1_ps(a31), m_b31=_mm_set1_ps(b31), m_c31=_mm_set1_ps(-c31);
			int k2=0;
			for(;k2+4<rgbn;k2+=4)
			{
				__m128 xb=_mm_loadu_ps(xbuffer+k2), yb=_mm_loadu_ps(ybuffer+k2), zb=_mm_loadu_ps(zbuffer+k2);
				__m128 den=_mm_mul_ps(m_a, xb);
				den=_mm_add_ps(den, _mm_mul_ps(m_b, yb));
				den=_mm_add_ps(den, _mm_mul_ps(m_c, zb));
				__m128 Xtx=_mm_mul_ps(m_B1, xb);
				__m128 Ytx=_mm_mul_ps(m_B4, xb);//y
				Xtx=_mm_add_ps(Xtx, _mm_mul_ps(m_B2, yb));
				Ytx=_mm_add_ps(Ytx, _mm_mul_ps(m_B5, yb));//y
				Xtx=_mm_add_ps(Xtx, _mm_mul_ps(m_B3, zb));
				Ytx=_mm_add_ps(Ytx, _mm_mul_ps(m_B6, zb));//y
				Xtx=_mm_div_ps(Xtx, den);
				Ytx=_mm_div_ps(Ytx, den);//y

				__m128 t1=_mm_mul_ps(den, m_t);
				__m128i cmp_t1=_mm_castps_si128(_mm_cmpge_ps(t1, _mm_setzero_ps()));

				__m128 l1=_mm_mul_ps(m_a12, Xtx);
				l1=_mm_add_ps(l1, _mm_mul_ps(m_b12, Ytx));
				__m128i cmp_l1=_mm_castps_si128(_mm_cmpgt_ps(l1, m_c12));

				__m128 l2=_mm_mul_ps(m_a23, Xtx);
				l2=_mm_add_ps(l2, _mm_mul_ps(m_b23, Ytx));
				__m128i cmp_l2=_mm_castps_si128(_mm_cmpgt_ps(l2, m_c23));

				__m128 l3=_mm_mul_ps(m_a31, Xtx);
				l3=_mm_add_ps(l3, _mm_mul_ps(m_b31, Ytx));
				__m128i cmp_l3=_mm_castps_si128(_mm_cmpgt_ps(l3, m_c31));

				__m128i cmp=_mm_and_si128(cmp_t1, cmp_l1);
				cmp=_mm_and_si128(cmp, cmp_l2);
				cmp=_mm_and_si128(cmp, cmp_l3);

				__m128i cmp2=_mm_shuffle_epi32(cmp, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
				__m128i cmp_or=_mm_or_si128(cmp, cmp2);//r3|r1, r2|r0, r1|r3, r0|r2
				cmp2=_mm_shuffle_epi32(cmp_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
				cmp_or=_mm_or_si128(cmp_or, cmp2);//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
				if(cmp_or.m128i_i32[0])
				{
					__m128 wb=_mm_loadu_ps(wbuffer+k2);
					__m128 m_adm=_mm_div_ps(den, m_t);
					__m128i c_a=_mm_castps_si128(_mm_cmpgt_ps(m_adm, wb));

					__m128i c_a2=_mm_shuffle_epi32(c_a, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
					__m128i c_a_or=_mm_or_si128(c_a, c_a2);//r3|r1, r2|r0, r1|r3, r0|r2
					c_a2=_mm_shuffle_epi32(c_a_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
					c_a_or=_mm_or_si128(c_a_or, c_a2);//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
					if(c_a_or.m128i_i32[0])
					{
						__m128i condition=_mm_and_si128(cmp, c_a);
						__m128i condition2=_mm_xor_si128(condition, _mm_set1_epi32(-1));

						wb=_mm_castsi128_ps(_mm_and_si128(condition2, _mm_castps_si128(wb)));
						__m128 wb2=_mm_castsi128_ps(_mm_and_si128(condition, _mm_castps_si128(m_adm)));
						wb2=_mm_castsi128_ps(_mm_or_si128(_mm_castps_si128(wb2), _mm_castps_si128(wb)));
						_mm_storeu_ps(wbuffer+k2, wb2);

						//mod operation		x%q = x-q*floor(x/q)
						__m128 Xtx_w=_mm_div_ps(Xtx, s_txwk);
						__m128 Ytx_h=_mm_div_ps(Ytx, s_txhk);

						//Xtx_w=_mm_floor_ps(Xtx_w);//floor(x/txw)	SSE4.1
						//Ytx_h=_mm_floor_ps(Ytx_h);
						Xtx_w=_mm_sub_ps(Xtx_w, half);//floor(x/q)
						Ytx_h=_mm_sub_ps(Ytx_h, half);
						__m128i i_Xtx_w=_mm_cvtps_epi32(Xtx_w);
						__m128i i_Ytx_h=_mm_cvtps_epi32(Ytx_h);
						Xtx_w=_mm_cvtepi32_ps(i_Xtx_w);
						Ytx_h=_mm_cvtepi32_ps(i_Ytx_h);

						Xtx_w=_mm_mul_ps(Xtx_w, s_txwk);//floor(x/txw)*txw
						Ytx_h=_mm_mul_ps(Ytx_h, s_txhk);

						Xtx=_mm_sub_ps(Xtx, Xtx_w);//x-floor(x/txw)*txw
						Ytx=_mm_sub_ps(Ytx, Ytx_h);
						
						//Xtx=_mm_floor_ps(Xtx);//99.9 -> 99 not 100	SSE4.1
						//Ytx=_mm_floor_ps(Ytx);
						Xtx=_mm_sub_ps(Xtx, half);
						Ytx=_mm_sub_ps(Ytx, half);

						__m128i i_Xtx=_mm_cvtps_epi32(Xtx);
						__m128i i_Ytx=_mm_cvtps_epi32(Ytx);

						//__m128i tx_idx=_mm_mullo_epi32(i_txwk, i_Ytx);//SSE4.1
						__m128i tx_idx_lo=_mm_mul_epu32(i_Ytx, i_txwk);
						i_Ytx=_mm_shuffle_epi32(i_Ytx, _MM_SHUFFLE(2, 3, 0, 1));//y2 y3 y0 y1
						__m128i tx_idx_hi=_mm_mul_epu32(i_Ytx, i_txwk);
						tx_idx_hi=_mm_shuffle_epi32(tx_idx_hi, _MM_SHUFFLE(2, 3, 0, 1));
						tx_idx_lo=_mm_and_si128(tx_idx_lo, mask2_lo);
						tx_idx_hi=_mm_and_si128(tx_idx_hi, mask2_hi);
						__m128i tx_idx=_mm_or_si128(tx_idx_lo, tx_idx_hi);

						tx_idx=_mm_add_epi32(tx_idx, i_Xtx);

						__m128i m_tx=_mm_set_epi32(txk[tx_idx.m128i_i32[3]], txk[tx_idx.m128i_i32[2]], txk[tx_idx.m128i_i32[1]], txk[tx_idx.m128i_i32[0]]);
						__m128i m_rgb=_mm_loadu_si128((__m128i*)(rgb+k2));
						__m128i m_rgb2=_mm_and_si128(m_rgb, condition2);
						m_tx=_mm_and_si128(m_tx, condition);
						m_rgb=_mm_and_si128(m_rgb, condition);

						__m128i const zero=_mm_setzero_si128();
						__m128i m_tx_lo=_mm_unpacklo_epi8(m_tx, zero);//{0, tx[7], ..., 0, tx[0]}
						__m128i m_tx_hi=_mm_unpackhi_epi8(m_tx, zero);//{0, tx[15], ..., 0, tx[8]}
						__m128i m_rgb_lo=_mm_unpacklo_epi8(m_rgb, zero);
						__m128i m_rgb_hi=_mm_unpackhi_epi8(m_rgb, zero);
						m_tx_lo=_mm_mullo_epi16(m_tx_lo, m_rgb_lo);//{{0, tx[7]}*{0, rgb[7]}, ..., {0, tx[0]}*{0, rgb[0]}}
						m_tx_hi=_mm_mullo_epi16(m_tx_hi, m_rgb_hi);
						m_tx_lo=_mm_shuffle_epi8(m_tx_lo, shuffle_lo);
						m_tx_hi=_mm_shuffle_epi8(m_tx_hi, shuffle_hi);
						m_tx_lo=_mm_and_si128(m_tx_lo, mask_lo);
						m_tx_hi=_mm_and_si128(m_tx_hi, mask_hi);
						m_tx=_mm_or_si128(m_tx_lo, m_tx_hi);

						m_tx=_mm_or_si128(m_tx, m_rgb2);
						_mm_storeu_si128((__m128i*)(rgb+k2), m_tx);
					}
				}
			}
			for(;k2<rgbn;++k2)
			{
				cpt=a*xbuffer[k2]+b*ybuffer[k2]+c*zbuffer[k2];
				float Xtx=(B1*xbuffer[k2]+B2*ybuffer[k2]+B3*zbuffer[k2])/cpt;
				float Ytx=(B4*xbuffer[k2]+B5*ybuffer[k2]+B6*zbuffer[k2])/cpt;
				if(cpt*t>=0&&a12*Xtx+b12*Ytx+c12>0&&a23*Xtx+b23*Ytx+c23>0&&a31*Xtx+b31*Ytx+c31>0)
				{
					admittance=cpt/t;
					if(admittance>wbuffer[k2])
					{
						int Xtx2=int(Xtx)%txwk;if(Xtx2<0)Xtx2+=txwk;
						int Ytx2=int(Ytx)%txhk;if(Ytx2<0)Ytx2+=txhk;
						((unsigned char*)&rgb[k2])[0]=((unsigned char*)&rgb[k2])[0]*((unsigned char*)&txk[Xtx2+Ytx2*txwk])[0]>>8;//b	//f(x,y)=x*y
						((unsigned char*)&rgb[k2])[1]=((unsigned char*)&rgb[k2])[1]*((unsigned char*)&txk[Xtx2+Ytx2*txwk])[1]>>8;//g
						((unsigned char*)&rgb[k2])[2]=((unsigned char*)&rgb[k2])[2]*((unsigned char*)&txk[Xtx2+Ytx2*txwk])[2]>>8;//r
					}
				}
			}
			//for(int k2=0;k2<rgbn;++k2)
			//{
			////	cpt=a*xbuffer[k2]+b*ybuffer[k2]+c*zbuffer[k2];//
			//	//auto &ux=xbuffer[k2], &uy=ybuffer[k2], &uz=zbuffer[k2];
			//	//cpt=a*ux+b*uy+c*zbuffer[k2];
			//	//float Xtx=(B1*ux+B2*uy+B3*uz)/cpt;
			//	//float Ytx=(B4*ux+B5*uy+B6*uz)/cpt;
			//	cpt=a*xbuffer[k2]+b*ybuffer[k2]+c*zbuffer[k2];
			//	float Xtx=(B1*xbuffer[k2]+B2*ybuffer[k2]+B3*zbuffer[k2])/cpt;
			//	float Ytx=(B4*xbuffer[k2]+B5*ybuffer[k2]+B6*zbuffer[k2])/cpt;
			//	if(cpt*t>=0&&a12*Xtx+b12*Ytx+c12>0&&a23*Xtx+b23*Ytx+c23>0&&a31*Xtx+b31*Ytx+c31>0)
			//	{
			//		admittance=cpt/t;
			//		if(admittance>wbuffer[k2])
			//		{
			//			int Xtx2=int(Xtx)%txwk;if(Xtx2<0)Xtx2+=txwk;
			//			int Ytx2=int(Ytx)%txhk;if(Ytx2<0)Ytx2+=txhk;
			//			((unsigned char*)&rgb[k2])[0]=((unsigned char*)&rgb[k2])[0]*((unsigned char*)&txk[Xtx2+Ytx2*txwk])[0]>>8;//b	//f(x,y)=x*y
			//			((unsigned char*)&rgb[k2])[1]=((unsigned char*)&rgb[k2])[1]*((unsigned char*)&txk[Xtx2+Ytx2*txwk])[1]>>8;//g
			//			((unsigned char*)&rgb[k2])[2]=((unsigned char*)&rgb[k2])[2]*((unsigned char*)&txk[Xtx2+Ytx2*txwk])[2]>>8;//r
			//		}
			//	}
			////	admittance=cpt/t;
			////	if(admittance>=0&&a12*Xtx+b12*Ytx+c12>0&&a23*Xtx+b23*Ytx+c23>0&&a31*Xtx+b31*Ytx+c31>0&&admittance>wbuffer[k2])
			//////	if(admittance>wbuffer[k2]&&a12*Xtx+b12*Ytx+c12>0&&a23*Xtx+b23*Ytx+c23>0&&a31*Xtx+b31*Ytx+c31>0)
			////	{
			////		int Xtx2=int(Xtx)%txwk;if(Xtx2<0)Xtx2+=txwk;
			////		int Ytx2=int(Ytx)%txhk;if(Ytx2<0)Ytx2+=txhk;
			////		((unsigned char*)&rgb[k2])[0]=((unsigned char*)&rgb[k2])[0]*((unsigned char*)&txk[Xtx2+Ytx2*txwk])[0]>>8;//b	//f(x,y)=x*y
			////		((unsigned char*)&rgb[k2])[1]=((unsigned char*)&rgb[k2])[1]*((unsigned char*)&txk[Xtx2+Ytx2*txwk])[1]>>8;//g
			////		((unsigned char*)&rgb[k2])[2]=((unsigned char*)&rgb[k2])[2]*((unsigned char*)&txk[Xtx2+Ytx2*txwk])[2]>>8;//r
			////	}
			//	//if(t*cpt>=0&&a12*Xtx+b12*Ytx+c12>0&&a23*Xtx+b23*Ytx+c23>0&&a31*Xtx+b31*Ytx+c31>0)
			//	//{
			//	//	admittance=cpt*cpt/t2;
			//	//	if(admittance>wbuffer[k2])
			//	//	{
			//	//		int Xtx2=int(Xtx)%txwk;if(Xtx2<0)Xtx2+=txwk;
			//	//		int Ytx2=int(Ytx)%txhk;if(Ytx2<0)Ytx2+=txhk;
			//	//	//	wbuffer[k2]=admittance;
			//	//	//	rgb[k2]=txk[Xtx2+Ytx2*txwk];
			//	//		((unsigned char*)&rgb[k2])[0]=((unsigned char*)&rgb[k2])[0]*((unsigned char*)&txk[Xtx2+Ytx2*txwk])[0]>>8;//b	//f(x,y)=x*y
			//	//		((unsigned char*)&rgb[k2])[1]=((unsigned char*)&rgb[k2])[1]*((unsigned char*)&txk[Xtx2+Ytx2*txwk])[1]>>8;//g
			//	//		((unsigned char*)&rgb[k2])[2]=((unsigned char*)&rgb[k2])[2]*((unsigned char*)&txk[Xtx2+Ytx2*txwk])[2]>>8;//r
			//	//	}
			//	//}
			//}
		}
	}
	else//transparent solid color
	{
		float a12=m12,											r12=a12*Xp3;		a12/=r12;
		float a23=Yp3-m12,	b23=-Xp3,	c23=-b23*m12,			r23=c23;			a23/=r23, b23/=r23, c23/=r23;
		float a31=-Yp3,		b31=Xp3,	c31=-a31*Xp3-b31*Yp3,	r31=b31*m12+c31;	a31/=r31, b31/=r31, c31/=r31;
		int color=tc;
	//	b/=a, c/=a, A1/=a, A2/=a, A3/=a, A4/=a, A5/=a, A6/=a;//cpt
		if(SSE4_1)
		{
			int k2=0;
			__m128i const shuffle_lo=_mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 15, 13, 11, 9, 7, 5, 3, 1);
			__m128i const shuffle_hi=_mm_set_epi8(15, 13, 11, 9, 7, 5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0);
			__m128i const mask_lo=_mm_set_epi32(0, 0, -1, -1);
			__m128i const mask_hi=_mm_set_epi32(-1, -1, 0, 0);
			__m128i m_color=_mm_set1_epi32(color);
			__m128 m_a=_mm_set1_ps(a), m_b=_mm_set1_ps(b), m_c=_mm_set1_ps(c);
			__m128 m_t=_mm_set1_ps(t);
			__m128 m_A1=_mm_set1_ps(A1), m_A2=_mm_set1_ps(A2), m_A3=_mm_set1_ps(A3);
			__m128 m_A4=_mm_set1_ps(A4), m_A5=_mm_set1_ps(A5), m_A6=_mm_set1_ps(A6);
			__m128 m_a12=_mm_set1_ps(a12);
			__m128 m_a23=_mm_set1_ps(a23), m_b23=_mm_set1_ps(b23), m_c23=_mm_set1_ps(-c23);
			__m128 m_a31=_mm_set1_ps(a31), m_b31=_mm_set1_ps(b31), m_c31=_mm_set1_ps(-c31);
			for(;k2+4<rgbn;k2+=4)
			{
				__m128 xb=_mm_loadu_ps(xbuffer+k2), yb=_mm_loadu_ps(ybuffer+k2), zb=_mm_loadu_ps(zbuffer+k2);
				__m128 den=_mm_mul_ps(m_a, xb);
				den=_mm_add_ps(den, _mm_mul_ps(m_b, yb));
				den=_mm_add_ps(den, _mm_mul_ps(m_c, zb));
				__m128 Xtx=_mm_mul_ps(m_A1, xb);
				__m128 Ytx=_mm_mul_ps(m_A4, xb);//y
				Xtx=_mm_add_ps(Xtx, _mm_mul_ps(m_A2, yb));
				Ytx=_mm_add_ps(Ytx, _mm_mul_ps(m_A5, yb));//y
				Xtx=_mm_add_ps(Xtx, _mm_mul_ps(m_A3, zb));
				Ytx=_mm_add_ps(Ytx, _mm_mul_ps(m_A6, zb));//y
				Xtx=_mm_div_ps(Xtx, den);
				Ytx=_mm_div_ps(Ytx, den);//y

				__m128 t1=_mm_mul_ps(den, m_t);
				__m128i cmp_t1=_mm_castps_si128(_mm_cmpge_ps(t1, _mm_setzero_ps()));

				__m128 l1=_mm_mul_ps(m_a12, Xtx);
				__m128i cmp_l1=_mm_castps_si128(_mm_cmpgt_ps(l1, _mm_setzero_ps()));

				__m128 l2=_mm_mul_ps(m_a23, Xtx);
				l2=_mm_add_ps(l2, _mm_mul_ps(m_b23, Ytx));
				__m128i cmp_l2=_mm_castps_si128(_mm_cmpgt_ps(l2, m_c23));

				__m128 l3=_mm_mul_ps(m_a31, Xtx);
				l3=_mm_add_ps(l3, _mm_mul_ps(m_b31, Ytx));
				__m128i cmp_l3=_mm_castps_si128(_mm_cmpgt_ps(l3, m_c31));

				__m128i cmp=_mm_and_si128(cmp_t1, cmp_l1);
				cmp=_mm_and_si128(cmp, cmp_l2);
				cmp=_mm_and_si128(cmp, cmp_l3);

				__m128i cmp2=_mm_shuffle_epi32(cmp, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
				__m128i cmp_or=_mm_or_si128(cmp, cmp2);//r3|r1, r2|r0, r1|r3, r0|r2
				cmp2=_mm_shuffle_epi32(cmp_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
				cmp_or=_mm_or_si128(cmp_or, cmp2);//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
				if(cmp_or.m128i_i32[0])
				{
					__m128 wb=_mm_loadu_ps(wbuffer+k2);
					__m128 m_adm=_mm_div_ps(den, m_t);
					__m128i c_a=_mm_castps_si128(_mm_cmpgt_ps(m_adm, wb));
					
					__m128i c_a2=_mm_shuffle_epi32(c_a, _MM_SHUFFLE(1, 0, 3, 2));//r1, r0, r3, r2
					__m128i c_a_or=_mm_or_si128(c_a, c_a2);//r3|r1, r2|r0, r1|r3, r0|r2
					c_a2=_mm_shuffle_epi32(c_a_or, _MM_SHUFFLE(0, 3, 2, 1));//r0|r2, r3|r1, r2|r0, r1|r3
					c_a_or=_mm_or_si128(c_a_or, c_a2);//r3|r1|r0|r2, r2|r0|r3|r1, r1|r3|r2|r0, r0|r2|r1|r3
					if(c_a_or.m128i_i32[0])
					{
						__m128i condition=_mm_and_si128(cmp, c_a);
						__m128i condition2=_mm_xor_si128(condition, _mm_set1_epi32(-1));

						wb=_mm_castsi128_ps(_mm_and_si128(condition2, _mm_castps_si128(wb)));
						__m128 wb2=_mm_castsi128_ps(_mm_and_si128(condition, _mm_castps_si128(m_adm)));
						wb2=_mm_castsi128_ps(_mm_or_si128(_mm_castps_si128(wb2), _mm_castps_si128(wb)));
						_mm_storeu_ps(wbuffer+k2, wb2);
						
						__m128i m_rgb=_mm_loadu_si128((__m128i*)(rgb+k2));
						__m128i m_color2=_mm_and_si128(m_color, condition);
						__m128i m_rgb2=_mm_and_si128(m_rgb, condition2);

						__m128i const zero=_mm_setzero_si128();
						__m128i m_color2_lo=_mm_unpacklo_epi8(m_color2, zero);//{0, tx[7], ..., 0, tx[0]}
						__m128i m_color2_hi=_mm_unpackhi_epi8(m_color2, zero);//{0, tx[15], ..., 0, tx[8]}
						__m128i m_rgb_lo=_mm_unpacklo_epi8(m_rgb, zero);
						__m128i m_rgb_hi=_mm_unpackhi_epi8(m_rgb, zero);
						m_color2_lo=_mm_mullo_epi16(m_color2_lo, m_rgb_lo);//{{0, tx[7]}*{0, rgb[7]}, ..., {0, tx[0]}*{0, rgb[0]}}
						m_color2_hi=_mm_mullo_epi16(m_color2_hi, m_rgb_hi);
						m_color2_lo=_mm_shuffle_epi8(m_color2_lo, shuffle_lo);
						m_color2_hi=_mm_shuffle_epi8(m_color2_hi, shuffle_hi);
						m_color2_lo=_mm_and_si128(m_color2_lo, mask_lo);
						m_color2_hi=_mm_and_si128(m_color2_hi, mask_hi);
						m_color2=_mm_or_si128(m_color2_lo, m_color2_hi);

						m_color2=_mm_or_si128(m_color2, m_rgb2);
						_mm_storeu_si128((__m128i*)(rgb+k2), m_color2);
					}
				}
			}
			for(;k2<rgbn;++k2)
			{
				cpt=a*xbuffer[k2]+b*ybuffer[k2]+c*zbuffer[k2];
				float Xp=(A1*xbuffer[k2]+A2*ybuffer[k2]+A3*zbuffer[k2])/cpt;
				float Yp=(A4*xbuffer[k2]+A5*ybuffer[k2]+A6*zbuffer[k2])/cpt;
				if(cpt*t>=0&&a12*Xp>0&&a23*Xp+b23*Yp+c23>0&&a31*Xp+b31*Yp+c31>0)
				{
					admittance=cpt/t;
					if(admittance>wbuffer[k2])
					{
						((unsigned char*)&rgb[k2])[0]=((unsigned char*)&rgb[k2])[0]*((unsigned char*)&color)[0]>>8;//b	//f(x,y)=x*y
						((unsigned char*)&rgb[k2])[1]=((unsigned char*)&rgb[k2])[1]*((unsigned char*)&color)[1]>>8;//g
						((unsigned char*)&rgb[k2])[2]=((unsigned char*)&rgb[k2])[2]*((unsigned char*)&color)[2]>>8;//r
					}
				}
			}
		}
		else
		{
			for(int k2=0;k2<rgbn;++k2)
			{
			//	cpt=xbuffer[k2]+b*ybuffer[k2]+c*zbuffer[k2];//
			/*	auto &ux=xbuffer[k2], &uy=ybuffer[k2], &uz=zbuffer[k2];
				cpt=a*ux+b*uy+c*zbuffer[k2];
				float Xp=(A1*ux+A2*uy+A3*uz)/cpt;
				float Yp=(A4*ux+A5*uy+A6*uz)/cpt;//*/
				cpt=a*xbuffer[k2]+b*ybuffer[k2]+c*zbuffer[k2];
				float Xp=(A1*xbuffer[k2]+A2*ybuffer[k2]+A3*zbuffer[k2])/cpt;
				float Yp=(A4*xbuffer[k2]+A5*ybuffer[k2]+A6*zbuffer[k2])/cpt;//*/
#if 1
				if(cpt*t>=0&&a12*Xp>0&&a23*Xp+b23*Yp+c23>0&&a31*Xp+b31*Yp+c31>0)
				{
					admittance=cpt/t;
					if(admittance>wbuffer[k2])
					{
						((unsigned char*)&rgb[k2])[0]=((unsigned char*)&rgb[k2])[0]*((unsigned char*)&color)[0]>>8;//b	//f(x,y)=x*y
						((unsigned char*)&rgb[k2])[1]=((unsigned char*)&rgb[k2])[1]*((unsigned char*)&color)[1]>>8;//g
						((unsigned char*)&rgb[k2])[2]=((unsigned char*)&rgb[k2])[2]*((unsigned char*)&color)[2]>>8;//r
					}
				}
										//	admittance=cpt/t;
		//	if(admittance>=0&&a12*Xp>0&&a23*Xp+b23*Yp+c23>0&&a31*Xp+b31*Yp+c31>0&&admittance>wbuffer[k2])
		////	if(admittance>wbuffer[k2]&&a12*Xp>0&&a23*Xp+b23*Yp+c23>0&&a31*Xp+b31*Yp+c31>0)
		//	{
		//		((unsigned char*)&rgb[k2])[0]=((unsigned char*)&rgb[k2])[0]*((unsigned char*)&color)[0]>>8;//b	//f(x,y)=x*y
		//		((unsigned char*)&rgb[k2])[1]=((unsigned char*)&rgb[k2])[1]*((unsigned char*)&color)[1]>>8;//g
		//		((unsigned char*)&rgb[k2])[2]=((unsigned char*)&rgb[k2])[2]*((unsigned char*)&color)[2]>>8;//r
		//	}
#else
				if(t*cpt>=0&&a12*Xp>0&&a23*Xp+b23*Yp+c23>0&&a31*Xp+b31*Yp+c31>0)
				{
					admittance=cpt*cpt/t2;
					if(admittance>wbuffer[k2])
					{
					//	wbuffer[k2]=admittance;
					//	rgb[k2]=color;
						((unsigned char*)&rgb[k2])[0]=((unsigned char*)&rgb[k2])[0]*((unsigned char*)&color)[0]>>8;//b	//f(x,y)=x*y
						((unsigned char*)&rgb[k2])[1]=((unsigned char*)&rgb[k2])[1]*((unsigned char*)&color)[1]>>8;//g
						((unsigned char*)&rgb[k2])[2]=((unsigned char*)&rgb[k2])[2]*((unsigned char*)&color)[2]>>8;//r
					}
				}
#endif
			}
		}
	}
}
void				NonlinearSW::render()
{
	for(int k=0, kEnd=objects.size();k<kEnd;++k)
	{
		auto &A=objects[k];
		for(int kt=0, ktEnd=A.tr.size();kt<ktEnd;++kt)
		{
			auto &tr=A.tr[kt];
			render_opaque(A.w[tr.a], A.w[tr.b], A.w[tr.c], tr.textured, tr.tc, tr.tx1, tr.tx2, tr.tx3, tr.txm);
		}
	}
	for(int k=0;k<npols;++k)
	{
		auto &p=env[k];
		render_opaque(p.p1, p.p2, p.p3, p.tx_idx>=0, p.tx_idx>=0?p.tx_idx:p.color, p.tx1, p.tx2, p.tx3, p.txm);
	}
	for(int k=0;k<ntpols;++k)
	{
		auto &p=tenv[k];
		render_transparent(p.p1, p.p2, p.p3, p.tx_idx>=0, p.tx_idx>=0?p.tx_idx:p.color, p.tx1, p.tx2, p.tx3, p.txm);
	}
	QueryPerformanceFrequency(&li);freq=li.QuadPart;
	QueryPerformanceCounter(&li);
	linelen=sprintf_s(line, "fps=%.10f, T=%.10fms, dcam=%.10f, fov=%.10f", freq/float(li.QuadPart-nticks), 1000.*(li.QuadPart-nticks)/freq, dcam, fov*720/_2pi);
	nticks=li.QuadPart;
	if(linelen>0)TextOutA(ghMemDC, 0, h-16, line, linelen);
}
void				NonlinearSW::show(){BitBlt(ghDC, 0, 0, w, h, ghMemDC, 0, 0, SRCCOPY);}
void				NonlinearSW::rebuffer1()//F2 cylindrical panorama
{
	float Xn, Yn, *ux_y=xbuffer, *uy_y=ybuffer, *uz_y=zbuffer;
	for(int ys=0;ys<h;++ys)
	{
		Yn=(ys-Y0)*tanfov/X0;
		float th=sqrt(1+Yn*Yn), uy_yx=Yn/th;
		for(int xs=0;xs<w;++xs)
		{
			Xn=(float(xs)/X0-1)*fov;
			ux_y[xs]=sin(Xn)/th, uy_y[xs]=uy_yx, uz_y[xs]=cos(Xn)/th;
		}
		ux_y=ux_y+w, uy_y=uy_y+w, uz_y=uz_y+w;
	}
}
void				NonlinearSW::rebuffer2()//F3 sideways cylindrical panorama
{
	float Xn, Yn, *ux_y=xbuffer, *uy_y=ybuffer, *uz_y=zbuffer;
	for(int ys=0;ys<h;++ys)
	{
		Yn=(ys-Y0)*fov/X0;
		float sinYn=sin(Yn), cosYn=cos(Yn);
		for(int xs=0;xs<w;++xs)
		{
			Xn=(float(xs)/X0-1)*tanfov;
			float th=sqrt(1+Xn*Xn);
			ux_y[xs]=Xn/th, uy_y[xs]=sinYn/th, uz_y[xs]=cosYn/th;
		}
		ux_y=ux_y+w, uy_y=uy_y+w, uz_y=uz_y+w;
	}
}
void				NonlinearSW::rebuffer3()//F4 spherical
{
	float Xn, Yn, *ux_y=xbuffer, *uy_y=ybuffer, *uz_y=zbuffer;
	//if(SSE4_1)
	//{
		const __m128i minus_one=_mm_set1_epi32(-1);
		const __m128 zero=_mm_setzero_ps(), one=_mm_set1_ps(1);
		__m128 m_X0=_mm_set1_ps(float(X0)), m_fov=_mm_set1_ps(fov);
		for(int ys=0;ys<h;++ys)
		{
			Yn=(ys-Y0)*fov/X0;
			int xs=0;
			__m128 m_Yn=_mm_set1_ps(Yn);
			__m128 m_Yn2=_mm_mul_ps(m_Yn, m_Yn);
			for(;xs+4<w;xs+=4)
			{
				__m128 m_Xn=_mm_div_ps(_mm_set_ps(float(xs+3), float(xs+2), float(xs+1), float(xs)), m_X0);
				m_Xn=_mm_sub_ps(m_Xn, one);
				m_Xn=_mm_mul_ps(m_Xn, m_fov);//Xn=(xs/X0-1)*fov

				__m128 m_th=_mm_mul_ps(m_Xn, m_Xn);
				m_th=_mm_add_ps(m_th, m_Yn2);//th=Xn^2+Yn^2

				__m128i th_cmp=_mm_castps_si128(_mm_cmpeq_ps(m_th, zero));
				__m128 m_ux=_mm_castsi128_ps(_mm_xor_si128(th_cmp, minus_one));
				__m128 m_uy=m_ux;
				__m128 m_uz=_mm_castsi128_ps(_mm_and_si128(th_cmp, _mm_castps_si128(one)));
				m_uz=_mm_or_ps(m_uz, m_ux);

				m_th=_mm_sqrt_ps(m_th);
				__m128 m_sr=_mm_sin_ps(m_th);
				//m_sr.m128_f32[0]=sin(m_th.m128_f32[0]);//
				//m_sr.m128_f32[1]=sin(m_th.m128_f32[1]);
				//m_sr.m128_f32[2]=sin(m_th.m128_f32[2]);
				//m_sr.m128_f32[3]=sin(m_th.m128_f32[3]);
				m_sr=_mm_div_ps(m_sr, m_th);

				__m128 m_ux2=_mm_mul_ps(m_Xn, m_sr);
				__m128 m_uy2=_mm_mul_ps(m_Yn, m_sr);
				__m128 m_uz2=_mm_cos_ps(m_th);
				m_ux=_mm_and_ps(m_ux, m_ux2);
				m_uy=_mm_and_ps(m_uy, m_uy2);
				m_uz=_mm_and_ps(m_uz, m_uz2);
				_mm_storeu_ps(ux_y+xs, m_ux);
				_mm_storeu_ps(uy_y+xs, m_uy);
				_mm_storeu_ps(uz_y+xs, m_uz);
			}
			for(;xs<w;++xs)
			{
				Xn=(float(xs)/X0-1)*fov;
				float th=sqrt(Xn*Xn+Yn*Yn);
				if(th==0){ux_y[xs]=uy_y[xs]=0, uz_y[xs]=1;continue;}
				float sr=sin(th)/th;
				ux_y[xs]=Xn*sr, uy_y[xs]=Yn*sr, uz_y[xs]=cos(th);
			}
			ux_y=ux_y+w, uy_y=uy_y+w, uz_y=uz_y+w;
		}
	//}
	//else//no SSE4.1
	//{
	//	for(int ys=0;ys<h;++ys)
	//	{
	//		Yn=(ys-Y0)*fov/X0;
	//		for(int xs=0;xs<w;++xs)
	//		{
	//			Xn=(float(xs)/X0-1)*fov;
	//			float th=sqrt(Xn*Xn+Yn*Yn);
	//			if(th==0){ux_y[xs]=uy_y[xs]=0, uz_y[xs]=1;continue;}
	//			float sr=sin(th)/th;
	//			ux_y[xs]=Xn*sr, uy_y[xs]=Yn*sr, uz_y[xs]=cos(th);
	//		}
	//		ux_y=ux_y+w, uy_y=uy_y+w, uz_y=uz_y+w;
	//	}
	//}
}
void				NonlinearSW::rebuffer4()//F5 fuzzy
{
	float Xn, Yn, *ux_y=xbuffer, *uy_y=ybuffer, *uz_y=zbuffer;
	for(int ys=0;ys<h;++ys)
	{
		for(int xs=0;xs<w;++xs)
		{
			Yn=(ys-Y0)*tanfov/X0;
			Xn=(float(xs)/X0-1)*tanfov;
										//	Yn=exp(Yn);

										//	Yn=(0.0000005-0.000001*rand())*tanfov+exp(Yn);
										//	Xn+=(0.0000005-0.000001*rand())*tanfov;

											Yn+=(0.0000005f-0.000001f*rand())*tanfov;
											Xn+=(0.0000005f-0.000001f*rand())*tanfov;

										//	Xn+=(0.0000001-0.0000002*rand())*tanfov	+0.05*exp(10*Xn);
										//	Xn+=0.1*cos(0.00001*rand());
										//	Xn+=0.05*cos(0.1*xs);
			float mag=sqrt(1+Xn*Xn+Yn*Yn);
			ux_y[xs]=Xn/mag, uy_y[xs]=Yn/mag, uz_y[xs]=1/mag;
											uy_y[xs]+=0.01f*cos(80*Xn);

			//Yn=(ys-Y0)*tanfov/X0;
			//Xn=(float(xs)/X0-1)*fov;
			//								Yn+=(0.000001-0.000002*rand())*tanfov;
			//								Xn+=(0.000001-0.000002*rand())*tanfov;
			//float r=sqrt(1+Xn*Xn+Yn*Yn);
			//ux_y[xs]=Xn/r, uy_y[xs]=Yn/r, uz_y[xs]=1/r;
		}
		ux_y=ux_y+w, uy_y=uy_y+w, uz_y=uz_y+w;
	}
}
void				NonlinearSW::pushTexture(int *texture_)
{
	++ntx;
	texture=(int**)realloc(texture, ntx*sizeof(int*)), texture[ntx-1]=texture_;
}
void				NonlinearSW::enqueueTextures(){}
int*				NonlinearSW::popTexture()
{
	ntx--;
	return texture[ntx];
}
void				NonlinearSW::clearTextures()
{
	for(int k=0;k<ntx;++k)
		free(texture[k]);
	ntx=0;
}
void				NonlinearSW::finish()
{
	DeleteObject(SelectObject(ghMemDC, hBitmap)), DeleteDC(ghMemDC);
	for(int k=0;k<ntx;++k)free(texture[k]);
	free(texture), free(wbuffer), free(xbuffer), free(ybuffer), free(zbuffer);
}
struct			ParallelSW:private Mode
{
	int **texture;
	int *libuffer, *lfbuffer;
	float *wbuffer;
	int *rgb, rgbn;
	HDC__ *ghMemDC;
	HBITMAP__ *hBitmap;
	void initiate();
	void resize();
	void changeFov();
	void clear_screen();
	void print(int x, int y, const char *a, ...);
	void print(Vector3f const &p, const char *a, ...);
	void draw_point(Vector3f const &p, int color);
	void draw_line(float x1, float y1, float x2, float y2, int lineColor);
	void draw_line(Vector3f const &p1, Vector3f const &p2, int lineColor);
	void draw_triangle(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, int color, bool transparent);
	void draw_ground();
	void render_opaque		(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, bool textured, int tc, Vector2f const &tx1, Matrix22 const &txm);
	void render_transparent	(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, bool textured, int tc, Vector2f const &tx1, Matrix22 const &txm);
	void render();
	void show();
	void draft_start		(float *v1, float *v2);
	void draft				(float *v1, float *v2);
	void pushTexture(int *texture_);
	void enqueueTextures();
	int* popTexture();
	void clearTextures();
	void finish();
} psw;
void				ParallelSW::initiate()
{
	texture=0;
	rgbn=h*w;
	wbuffer=(float*)malloc(rgbn*sizeof(float)), libuffer=(int*)malloc(h*sizeof(int)), lfbuffer=(int*)malloc(h*sizeof(int));
	ghMemDC=CreateCompatibleDC(ghDC);
	tagBITMAPINFO bmpInfo={{sizeof(tagBITMAPINFOHEADER), w, -h, 1, 32, BI_RGB, 0, 0, 0, 0, 0}};
	hBitmap=(HBITMAP__*)SelectObject(ghMemDC, CreateDIBSection(0, &bmpInfo, DIB_RGB_COLORS, (void**)&rgb, 0, 0));
	SetBkMode(ghMemDC, TRANSPARENT);
}
void				ParallelSW::resize()
{
	rgbn=h*w;
	wbuffer=(float*)realloc(wbuffer, rgbn*sizeof(float)), libuffer=(int*)realloc(libuffer, h*sizeof(int)), lfbuffer=(int*)realloc(lfbuffer, h*sizeof(int));
	DeleteObject((HBITMAP__*)SelectObject(ghMemDC, hBitmap));
	tagBITMAPINFO bmpInfo={{sizeof(tagBITMAPINFOHEADER), w, -h, 1, 32, BI_RGB, 0, 0, 0, 0, 0}};
	hBitmap=(HBITMAP__*)SelectObject(ghMemDC, CreateDIBSection(0, &bmpInfo, DIB_RGB_COLORS, (void**)&rgb, 0, 0));
	SetBkMode(ghMemDC, TRANSPARENT);
}
void				ParallelSW::changeFov(){}
void				ParallelSW::clear_screen()
{
	memset(rgb, 0xFF, rgbn*sizeof(int));
	for(int k=0;k<rgbn;++k)
		wbuffer[k]=infinity;
}
void				ParallelSW::print(int x, int y, const char *a, ...){GUIVPrint(ghMemDC, x, y, a, (char*)(&a+1));}
void				ParallelSW::print(Vector3f const &p, const char *a, ...)
{
	Vector2f s;
	Vector3f c;
	world_to_screen(p, s, c);
	if(c.z>0)
		GUIVPrint(ghMemDC, (int)s.x, (int)s.y, a, (char*)(&a+1));
}
void				ParallelSW::draw_point(Vector3f const &p, int color){}
void				ParallelSW::draw_line(float x1, float y1, float x2, float y2, int lineColor){}
void				ParallelSW::draw_line(Vector3f const &p1, Vector3f const &p2, int lineColor){}
void				ParallelSW::draw_triangle(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, int color, bool transparent)
{
	Vector2f zero;
	if(transparent)
		render_transparent(p1, p2, p3, false, color, zero, Matrix22());
	else
		render_opaque(p1, p2, p3, false, color, zero, Matrix22());
}
void				ParallelSW::draw_ground()
{
}
void				ParallelSW::render_opaque		(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, bool textured, int tc, Vector2f const &tx1, Matrix22 const &txm)
{
	float dx, dy, dz, cpt, v1[5], v2[5], v3[5], admittance;
	dx=p1.x-cam.x, dy=p1.y-cam.y, dz=p1.z-cam.z, cpt=dx*cax+dy*sax, v1[2]=dx*sax-dy*cax, v1[3]=cpt*say-dz*cay, v1[4]=cpt*cay+dz*say, cpt=X0/(1000*tanfov), v1[0]=X0+v1[2]*cpt, v1[1]=Y0+v1[3]*cpt;
	dx=p2.x-cam.x, dy=p2.y-cam.y, dz=p2.z-cam.z, cpt=dx*cax+dy*sax, v2[2]=dx*sax-dy*cax, v2[3]=cpt*say-dz*cay, v2[4]=cpt*cay+dz*say, cpt=X0/(1000*tanfov), v2[0]=X0+v2[2]*cpt, v2[1]=Y0+v2[3]*cpt;
	dx=p3.x-cam.x, dy=p3.y-cam.y, dz=p3.z-cam.z, cpt=dx*cax+dy*sax, v3[2]=dx*sax-dy*cax, v3[3]=cpt*say-dz*cay, v3[4]=cpt*cay+dz*say, cpt=X0/(1000*tanfov), v3[0]=X0+v3[2]*cpt, v3[1]=Y0+v3[3]*cpt;
	if(v1[4]<0&&v2[4]<0&&v3[4]<0)		return;
	for(int k2=0;k2<h;++k2)libuffer[k2]=w; memset(lfbuffer, 0, h*sizeof(int));
	draft_start(v1, v2), draft(v2, v3), draft(v3, v1);
	float au12=v2[2]-v1[2],				bu12=v2[3]-v1[3],			cu12=v2[4]-v1[4];									//u12	=<12>
	float aux3=v3[2]-v1[2],				bux3=v3[3]-v1[3],			cux3=v3[4]-v1[4];									//ux3	=<13>
	float a=bu12*cux3-bux3*cu12,		b=aux3*cu12-au12*cux3,		c=au12*bux3-aux3*bu12;								//abc	=<n>	=<12>x<13>
	cpt=1000*tanfov;
	float A1=-cpt*a/c/X0, A2=-cpt*b/c/X0, A3=a/c*(cpt+v1[2])+b/c*(cpt*Y0/X0+v1[3])+v1[4];
	if(textured)
	{
		float	ux3_1=sqrt(au12*au12+bu12*bu12+cu12*cu12);	au12/=ux3_1,		bu12/=ux3_1,		cu12/=ux3_1;		//u12	=<u12>	=<12>/|12|
				ux3_1=aux3*au12+bux3*bu12+cux3*cu12,		aux3-=ux3_1*au12,	bux3-=ux3_1*bu12,	cux3-=ux3_1*cu12;	//ux3	=<x3>	=<13>-<u12>.<13><u12>
				ux3_1=sqrt(aux3*aux3+bux3*bux3+cux3*cux3),	aux3/=ux3_1,		bux3/=ux3_1,		cux3/=ux3_1;		//ux3	=<ux3>	=<x3>/|x3|
				ux3_1=aux3*v1[2]+bux3*v1[3]+cux3*v1[4];																	//ux3_1			=<ux3>.<1>
		float	u12_1=au12*v1[2]+bu12*v1[3]+cu12*v1[4];																	//u12_1			=<u12>.<1>
		float	C1=aux3*txm.a+au12*txm.b, C2=bux3*txm.a+bu12*txm.b, C3=cux3*txm.a+cu12*txm.b;
		float	C4=aux3*txm.c+au12*txm.d, C5=bux3*txm.c+bu12*txm.d, C6=cux3*txm.c+cu12*txm.d;
		float	D1=C1*cpt/X0+C3*A1, D2=C2*cpt/X0+C3*A2, D3=-C1*cpt-C2*cpt*Y0/X0+C3*A3+tx1.x-ux3_1*txm.a-u12_1*txm.b;
		float	D4=C4*cpt/X0+C6*A1, D5=C5*cpt/X0+C6*A2, D6=-C4*cpt-C5*cpt*Y0/X0+C6*A3+tx1.y-ux3_1*txm.c-u12_1*txm.d;
		float *wbk=wbuffer;
		int libk, lfbk, *rgbk=rgb, *txk=texture[tc], txhk=txh[tc], txwk=txw[tc];
		for(int k2=0;k2<h;++k2)
		{
			libk=libuffer[k2]<0?0:libuffer[k2], lfbk=lfbuffer[k2]>w?w:lfbuffer[k2];
			admittance=A1*libk+A2*k2+A3, a=D1*libk+D2*k2+D3, b=D4*libk+D5*k2+D6;
			for(int k3=libk;k3<lfbk;++k3)
			{
				if(admittance>0&&admittance<wbk[k3])
				{
					wbk[k3]=admittance;
					int Xtx=int(a)%txwk;if(Xtx<0)Xtx+=txwk;
					int Ytx=int(b)%txhk;if(Ytx<0)Ytx+=txhk;
					rgbk[k3]=txk[Xtx+Ytx*txwk];
				}
				admittance+=A1, a+=D1, b+=D4;
			}
			wbk=wbk+w, rgbk=rgbk+w;
		}
	}
	else
	{
		float *wbk=wbuffer;
		int libk, lfbk, *rgbk=rgb, color=tc;
		for(int k2=0;k2<h;++k2)
		{
			libk=libuffer[k2]<0?0:libuffer[k2], lfbk=lfbuffer[k2]>w?w:lfbuffer[k2];
			admittance=A1*libk+A2*k2+A3;
			for(int k3=libk;k3<lfbk;++k3)
			{
				if(admittance>0&&admittance<wbk[k3])
					wbk[k3]=admittance, rgbk[k3]=color;
				admittance+=A1;
			}
			wbk=wbk+w, rgbk=rgbk+w;
		}
	}
}
void				ParallelSW::render_transparent	(Vector3f const &p1, Vector3f const &p2, Vector3f const &p3, bool textured, int tc, Vector2f const &tx1, Matrix22 const &txm)
{
	float dx, dy, dz, cpt, v1[5], v2[5], v3[5], admittance;
	dx=p1.x-cam.x, dy=p1.y-cam.y, dz=p1.z-cam.z, cpt=dx*cax+dy*sax, v1[2]=dx*sax-dy*cax, v1[3]=cpt*say-dz*cay, v1[4]=cpt*cay+dz*say, cpt=X0/(1000*tanfov), v1[0]=X0+v1[2]*cpt, v1[1]=Y0+v1[3]*cpt;
	dx=p2.x-cam.x, dy=p2.y-cam.y, dz=p2.z-cam.z, cpt=dx*cax+dy*sax, v2[2]=dx*sax-dy*cax, v2[3]=cpt*say-dz*cay, v2[4]=cpt*cay+dz*say, cpt=X0/(1000*tanfov), v2[0]=X0+v2[2]*cpt, v2[1]=Y0+v2[3]*cpt;
	dx=p3.x-cam.x, dy=p3.y-cam.y, dz=p3.z-cam.z, cpt=dx*cax+dy*sax, v3[2]=dx*sax-dy*cax, v3[3]=cpt*say-dz*cay, v3[4]=cpt*cay+dz*say, cpt=X0/(1000*tanfov), v3[0]=X0+v3[2]*cpt, v3[1]=Y0+v3[3]*cpt;
	if(v1[4]<0&&v2[4]<0&&v3[4]<0)		return;
	for(int k2=0;k2<h;++k2)libuffer[k2]=w; memset(lfbuffer, 0, h*sizeof(int));
	draft_start(v1, v2), draft(v2, v3), draft(v3, v1);
	float au12=v2[2]-v1[2],				bu12=v2[3]-v1[3],			cu12=v2[4]-v1[4];									//u12	=<12>
	float aux3=v3[2]-v1[2],				bux3=v3[3]-v1[3],			cux3=v3[4]-v1[4];									//ux3	=<13>
	float a=bu12*cux3-bux3*cu12,		b=aux3*cu12-au12*cux3,		c=au12*bux3-aux3*bu12;								//abc	=<n>	=<12>x<13>
	cpt=1000*tanfov;
	float A1=-cpt*a/c/X0, A2=-cpt*b/c/X0, A3=a/c*(cpt+v1[2])+b/c*(cpt*Y0/X0+v1[3])+v1[4];
	if(textured)
	{
		float	ux3_1=sqrt(au12*au12+bu12*bu12+cu12*cu12);	au12/=ux3_1,		bu12/=ux3_1,		cu12/=ux3_1;		//u12	=<u12>	=<12>/|12|
				ux3_1=aux3*au12+bux3*bu12+cux3*cu12,		aux3-=ux3_1*au12,	bux3-=ux3_1*bu12,	cux3-=ux3_1*cu12;	//ux3	=<x3>	=<13>-<u12>.<13><u12>
				ux3_1=sqrt(aux3*aux3+bux3*bux3+cux3*cux3),	aux3/=ux3_1,		bux3/=ux3_1,		cux3/=ux3_1;		//ux3	=<ux3>	=<x3>/|x3|
				ux3_1=aux3*v1[2]+bux3*v1[3]+cux3*v1[4];																	//ux3_1			=<ux3>.<1>
		float	u12_1=au12*v1[2]+bu12*v1[3]+cu12*v1[4];																	//u12_1			=<u12>.<1>
		float	C1=aux3*txm.a+au12*txm.b, C2=bux3*txm.a+bu12*txm.b, C3=cux3*txm.a+cu12*txm.b;
		float	C4=aux3*txm.c+au12*txm.d, C5=bux3*txm.c+bu12*txm.d, C6=cux3*txm.c+cu12*txm.d;
		float	D1=C1*cpt/X0+C3*A1, D2=C2*cpt/X0+C3*A2, D3=-C1*cpt-C2*cpt*Y0/X0+C3*A3+tx1.x-ux3_1*txm.a-u12_1*txm.b;
		float	D4=C4*cpt/X0+C6*A1, D5=C5*cpt/X0+C6*A2, D6=-C4*cpt-C5*cpt*Y0/X0+C6*A3+tx1.y-ux3_1*txm.c-u12_1*txm.d;
		float *wbk=wbuffer;
		int libk, lfbk, *rgbk=rgb, *txk=texture[tc], txhk=txh[tc], txwk=txw[tc];
		for(int k2=0;k2<h;++k2)
		{
			libk=libuffer[k2]<0?0:libuffer[k2], lfbk=lfbuffer[k2]>w?w:lfbuffer[k2];
			admittance=A1*libk+A2*k2+A3, a=D1*libk+D2*k2+D3, b=D4*libk+D5*k2+D6;
			for(int k3=libk;k3<lfbk;++k3)
			{
				if(admittance>0&&admittance<wbk[k3])
				{
					wbk[k3]=admittance;
					int Xtx=int(a)%txwk;if(Xtx<0)Xtx+=txwk;
					int Ytx=int(b)%txhk;if(Ytx<0)Ytx+=txhk;
					rgbk[k3]=txk[Xtx+Ytx*txwk];
				}
				admittance+=A1, a+=D1, b+=D4;
			}
			wbk=wbk+w, rgbk=rgbk+w;
		}
	}
	else
	{
		float *wbk=wbuffer;
		int libk, lfbk, *rgbk=rgb, color=tc;
		for(int k2=0;k2<h;++k2)
		{
			libk=libuffer[k2]<0?0:libuffer[k2], lfbk=lfbuffer[k2]>w?w:lfbuffer[k2];
			admittance=A1*libk+A2*k2+A3;
			for(int k3=libk;k3<lfbk;++k3)
			{
				if(admittance>0&&admittance<wbk[k3])
					wbk[k3]=admittance, rgbk[k3]=color;
				admittance+=A1;
			}
			wbk=wbk+w, rgbk=rgbk+w;
		}
	}
}
void				ParallelSW::render()
{
	for(int k=0, kEnd=objects.size();k<kEnd;++k)
	{
		auto &A=objects[k];
		for(int kt=0, ktEnd=A.tr.size();kt<ktEnd;++kt)
		{
			auto &tr=A.tr[kt];
			render_opaque(A.w[tr.a], A.w[tr.b], A.w[tr.c], tr.textured, tr.tc, tr.tx1, tr.txm);
		}
	}
	for(int k=0;k<npols;++k)
	{
		auto &p=env[k];
		render_opaque(p.p1, p.p2, p.p3, p.tx_idx>=0, p.tx_idx>=0?p.tx_idx:p.color, p.tx1, p.txm);
	}
	for(int k=0;k<ntpols;++k)
	{
		auto &p=tenv[k];
		render_transparent(p.p1, p.p2, p.p3, p.tx_idx>=0, p.tx_idx>=0?p.tx_idx:p.color, p.tx1, p.txm);
	}
	QueryPerformanceFrequency(&li);freq=li.QuadPart;
	QueryPerformanceCounter(&li);
	linelen=sprintf_s(line, 128, "fps=%.10f, T=%.10fms, dcam=%.10f, fov=%.10f", freq/float(li.QuadPart-nticks), 1000.*(li.QuadPart-nticks)/freq, dcam, atan(tanfov)*720/_2pi);
	nticks=li.QuadPart;
	if(linelen>0)TextOutA(ghMemDC, 0, h-16, line, linelen);
}
void				ParallelSW::show(){BitBlt(ghDC, 0, 0, w, h, ghMemDC, 0, 0, SRCCOPY);}
void				ParallelSW::draft_start		(float *v1, float *v2)
{
	int k2;float k3, A=(v2[0]-v1[0])/(v2[1]-v1[1]);
		 if(v1[1]<v2[1]){	k3=A*((v1[1]<0?0:v1[1])-v1[1])+v1[0];		 if(v1[0]<v2[0])for(long long k=long long(v1[1])<0?	0:long long(v1[1])	;k<h&&k<v2[1]	;++k)	k2=int(k3), k2=k2<0?v1[0]>0?int(v1[0]):0:k2>w?v2[0]<w?int(v2[0]):w:k2<v1[0]?int(v1[0]):k2>v2[0]?int(v2[0]):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;
																	else				for(long long k=long long(v1[1])<0?	0:long long(v1[1])	;k<h&&k<v2[1]	;++k)	k2=int(k3), k2=k2<0?v2[0]>0?int(v2[0]):0:k2>w?v1[0]<w?int(v1[0]):w:k2<v2[0]?int(v2[0]):k2>v1[0]?int(v1[0]):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;}
	else				{	k3=A*((v2[1]<0?0:v2[1])-v1[1])+v1[0];		 if(v1[0]<v2[0])for(long long k=long long(v2[1])<0?	0:long long(v2[1])	;k<h&&k<v1[1]	;++k)	k2=int(k3), k2=k2<0?v1[0]>0?int(v1[0]):0:k2>w?v2[0]<w?int(v2[0]):w:k2<v1[0]?int(v1[0]):k2>v2[0]?int(v2[0]):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;
																	else				for(long long k=long long(v2[1])<0?	0:long long(v2[1])	;k<h&&k<v1[1]	;++k)	k2=int(k3), k2=k2<0?v2[0]>0?int(v2[0]):0:k2>w?v1[0]<w?int(v1[0]):w:k2<v2[0]?int(v2[0]):k2>v1[0]?int(v1[0]):k2, libuffer[k]=lfbuffer[k]=k2, k3+=A;}
}
void				ParallelSW::draft			(float *v1, float *v2)
{
	int k2;float k3, A=(v2[0]-v1[0])/(v2[1]-v1[1]);
		 if(v1[1]<v2[1]){	k3=A*((v1[1]<0?0:v1[1])-v1[1])+v1[0];		 if(v1[0]<v2[0])for(long long k=long long(v1[1])<0?	0:long long(v1[1])	;k<h&&k<v2[1]	;++k)	k2=int(k3), k2=k2<0?v1[0]>0?int(v1[0]):0:k2>w?v2[0]<w?int(v2[0]):w:k2<v1[0]?int(v1[0]):k2>v2[0]?int(v2[0]):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;
																	else				for(long long k=long long(v1[1])<0?	0:long long(v1[1])	;k<h&&k<v2[1]	;++k)	k2=int(k3), k2=k2<0?v2[0]>0?int(v2[0]):0:k2>w?v1[0]<w?int(v1[0]):w:k2<v2[0]?int(v2[0]):k2>v1[0]?int(v1[0]):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;}
	else				{	k3=A*((v2[1]<0?0:v2[1])-v1[1])+v1[0];		 if(v1[0]<v2[0])for(long long k=long long(v2[1])<0?	0:long long(v2[1])	;k<h&&k<v1[1]	;++k)	k2=int(k3), k2=k2<0?v1[0]>0?int(v1[0]):0:k2>w?v2[0]<w?int(v2[0]):w:k2<v1[0]?int(v1[0]):k2>v2[0]?int(v2[0]):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;
																	else				for(long long k=long long(v2[1])<0?	0:long long(v2[1])	;k<h&&k<v1[1]	;++k)	k2=int(k3), k2=k2<0?v2[0]>0?int(v2[0]):0:k2>w?v1[0]<w?int(v1[0]):w:k2<v2[0]?int(v2[0]):k2>v1[0]?int(v1[0]):k2, libuffer[k]=k2>=libuffer[k]?libuffer[k]:k2, lfbuffer[k]=k2<=lfbuffer[k]?lfbuffer[k]:k2, k3+=A;}
}
void				ParallelSW::pushTexture(int *texture_)
{
	++ntx;
	texture=(int**)realloc(texture, ntx*sizeof(int*)), texture[ntx-1]=texture_;
}
void				ParallelSW::enqueueTextures(){}
int*				ParallelSW::popTexture()
{
	ntx--;
	return texture[ntx];
}
void				ParallelSW::clearTextures()
{
	for(int k=0;k<ntx;++k)
		free(texture[k]);
	ntx=0;
}
void				ParallelSW::finish()
{
	DeleteObject(SelectObject(ghMemDC, hBitmap)), DeleteDC(ghMemDC);
	for(int k=0;k<ntx;++k)free(texture[k]);
	free(texture), free(libuffer), free(lfbuffer), free(wbuffer);
}
std::wstring	wider(char *a)
{
	auto len=strlen(a);
	auto buf=new wchar_t[len+1];
	swprintf_s(buf, len+1, L"%S", a);
	std::wstring str=buf;
	delete[] buf;
	return str;
}
std::string		narrower(std::wstring str)
{
	auto len=str.size();
	auto buf=new char[len+1];
	sprintf_s(buf, len+1, "%S", str.c_str());
	std::string str2=buf;
	delete[] buf;
	return str2;
}
bool			loadText(std::string &str)
{
	struct stat info;
	if(!stat(str.c_str(), &info))
	{
		std::ifstream file(str.c_str(), std::ios::in);
		if(file.is_open())
		{
			char *a=(char*)malloc((info.st_size+1)*sizeof(char));
			file.read(a, info.st_size), a[info.st_size]='\0';
			file.close();
			str=a;
			return true;
		}
	}
	str="";
	return false;
}
int				comtok(std::string &str, int *pI, int *pF)
{
	for(*pI=*pF;*pI==*pF;)
	{
		if(*pF<long long(str.size()))
		{
			if(str[*pF]=='/')
			{
				if(*pF+1<long long(str.size())&&str[*pF+1]=='*')
				{
					for(int k=*pF+2;;++k)
					{
						if(k==long long(str.size()))
						{
							*pI=k;
							break;
						}
						if(str[k]=='*'&&k+1<long long(str.size())&&str[k+1]=='/')
						{
							*pI=k+2;
							break;
						}
					}
				}
				else
				{
					for(int k=*pF+1;;++k)
					{
						if(k==long long(str.size()))
						{
							*pI=k;
							break;
						}
						if(str[k]=='\r'&&k+1<long long(str.size())&&str[k+1]=='\n')
						{
							*pI=k+2;
							break;
						}
						else if(str[k]=='\r'||str[k]=='\n')
						{
							*pI=k+1;
							break;
						}
					}
				}
			}
			else if(str[*pF]=='\r'&&*pF+1<long long(str.size())&&str[*pF+1]=='\n')
			{
				*pI=*pF+2;
			}
			else if(str[*pF]=='\r'||str[*pF]=='\n')
			{
				*pI=*pF+1;
			}
		}
		else
		{
			*pI=*pF;
		}
		if(*pI==long long(str.size()))
			return 0;
		for(int k=*pI;;++k)
		{
			if(k==long long(str.size())||str[k]=='\r'||str[k]=='\n'||str[k]=='/')
			{
				*pF=k;
				break;
			}
		}
	}
	return 1;
}
float			readDoubleFromBuffer(std::string &str, int i, int f, int *pF=0, int *pPrev=0)
{
	if(pF)*pF=-1;
	if(pPrev)*pPrev=0;
	int n=1;
	float value=0;
	if(str.size())
	{
		for(int k=i;k<f;++k)
		{
			switch(str[k])
			{
			case '-':
				n*=-1;
				continue;
			case '1':case '2':case '3':case '4':case '5':case '6':case '7':case '8':case '9':case '0':case '.':
				{
					int k2=k+1, endValid=0;
					for(;k2<f;++k2)
						if((str[k2]<'0'||str[k2]>'9')&&str[k2]!='.')
						{
							endValid=1;
							break;
						}
					if(!endValid)
						k2=f;
					if(pF)*pF=k2==f?-1:k2;
					for(int k3=k;k3<k2;++k3)
					{
						if(str[k3]>='0'&&str[k3]<='9')
						{
							float p=1;
							for(int k4=k;k4<k2;++k4)
							{
								if(str[k4]=='.')
								{
									for(int k5=k4+1;k5<k2;++k5)
									{
										if(str[k5]=='.')continue;
										p/=10;
									}
									break;
								}
							}
							for(int k4=k2-1;k4>=k;k4--)
							{
								if(str[k4]=='.')continue;
								value+=(str[k4]-'0')*p, p*=10;
							}
							break;
						}
					}
				}
				break;
			case '`':
				if(pPrev)*pPrev=1;
				if(pF)*pF=k+1;
				break;
			default:
				continue;
			}
			break;
		}
	}
	return n*value;
}
void			transform(Triangle0 &p)
{
//	auto &p=env[k];
	Vector3f _12=p.p2-p.p1, _13=p.p3-p.p1;
	float _r=_12.dot(_13)/_12.mag_sq();
	Vector3f _x3=_13-_r*_12;
	float im_x3=1/_x3.magnitude(), im_12=1/_12.magnitude();
	Vector2f tx12=p.tx2-p.tx1;
	p.txm.a=(p.tx3.x-p.tx1.x-_r*tx12.x)*im_x3, p.txm.b=tx12.x*im_12;
	p.txm.c=(p.tx3.y-p.tx1.y-_r*tx12.y)*im_x3, p.txm.d=tx12.y*im_12;
	//float _x12=p.p2.x-p.p1.x, _y12=p.p2.y-p.p1.y, _z12=p.p2.z-p.p1.z;
	//float _x13=p.p3.x-p.p1.x, _y13=p.p3.y-p.p1.y, _z13=p.p3.z-p.p1.z;
	//float _r=(_x12*_x13+_y12*_y13+_z12*_z13)/(_x12*_x12+_y12*_y12+_z12*_z12);
	//float _xx3=_x13-_r*_x12, _yx3=_y13-_r*_y12, _zx3=_z13-_r*_z12;
	//float _x3=sqrt(_xx3*_xx3+_yx3*_yx3+_zx3*_zx3), _12=sqrt(_x12*_x12+_y12*_y12+_z12*_z12);
	//p.txm.a=(p.tx3.x-p.tx1.x-_r*(p.tx2.x-p.tx1.x))/_x3, p.txm.b=(p.tx2.x-p.tx1.x)/_12;
	//p.txm.c=(p.tx3.y-p.tx1.y-_r*(p.tx2.y-p.tx1.y))/_x3, p.txm.d=(p.tx2.y-p.tx1.y)/_12;
}
//void			transform_transparent(int k)
//{
//	auto &p=tenv[k];
//	Vector3f _12=p.p2-p.p1, _13=p.p3-p.p1;
//	float _r=_12.dot(_13)/_12.mag_sq();
//	Vector3f _x3=_13-_r*_12;
//	float im_x3=1/_x3.magnitude(), im_12=1/_12.magnitude();
//	Vector2f tx12=p.tx2-p.tx1;
//	p.txm.a=(p.tx3.x-p.tx1.x-_r*tx12.x)*im_x3, p.txm.b=tx12.x*im_12;
//	p.txm.c=(p.tx3.y-p.tx1.y-_r*tx12.y)*im_x3, p.txm.d=tx12.y*im_12;
//	//float _x12=p.p2.x-p.p1.x, _y12=p.p2.y-p.p1.y, _z12=p.p2.z-p.p1.z;
//	//float _x13=p.p3.x-p.p1.x, _y13=p.p3.y-p.p1.y, _z13=p.p3.z-p.p1.z;
//	//float _r=(_x12*_x13+_y12*_y13+_z12*_z13)/(_x12*_x12+_y12*_y12+_z12*_z12);
//	//float _xx3=_x13-_r*_x12, _yx3=_y13-_r*_y12, _zx3=_z13-_r*_z12;
//	//float _x3=sqrt(_xx3*_xx3+_yx3*_yx3+_zx3*_zx3), _12=sqrt(_x12*_x12+_y12*_y12+_z12*_z12);
//	//p.txm.a=(p.tx3.x-p.tx1.x-_r*(p.tx2.x-p.tx1.x))/_x3, p.txm.b=(p.tx2.x-p.tx1.x)/_12;
//	//p.txm.c=(p.tx3.y-p.tx1.y-_r*(p.tx2.y-p.tx1.y))/_x3, p.txm.d=(p.tx2.y-p.tx1.y)/_12;
//}
void			antitransform(int k)
{
	auto &p=env[k];
	float _x12=p.p2.x-p.p1.x, _y12=p.p2.y-p.p1.y, _z12=p.p2.z-p.p1.z;
	float _x13=p.p3.x-p.p1.x, _y13=p.p3.y-p.p1.y, _z13=p.p3.z-p.p1.z;
	float _r=(_x12*_x13+_y12*_y13+_z12*_z13)/(_x12*_x12+_y12*_y12+_z12*_z12);
	float _xx3=_x13-_r*_x12, _yx3=_y13-_r*_y12, _zx3=_z13-_r*_z12;
	float _x3=sqrt(_xx3*_xx3+_yx3*_yx3+_zx3*_zx3), _12=sqrt(_x12*_x12+_y12*_y12+_z12*_z12);

	float	p2x=0,		p2y=_12,
			p3x=_x3,	p3y=(_x12*_x13+_y12*_y13+_z12*_z13)/_12;
	p.tx2.x=p.txm.a*p2x+p.txm.b*p2y, p.tx2.y=p.txm.c*p2x+p.txm.d*p2y;
	p.tx3.x=p.txm.a*p3x+p.txm.b*p3y, p.tx3.y=p.txm.c*p3x+p.txm.d*p3y;
//	p.tx1.x=		p.tx1.y=0;
//	p.tx2.x=0,	p.tx2.y=_12;
//	p.tx3.x=_x3,	p.tx3.y=(_x12*_x13+_y12*_y13+_z12*_z13)/_12;
}
void			antitransform_transparent(int k)
{
	auto &p=tenv[k];
	float _x12=p.p2.x-p.p1.x, _y12=p.p2.y-p.p1.y, _z12=p.p2.z-p.p1.z;
	float _x13=p.p3.x-p.p1.x, _y13=p.p3.y-p.p1.y, _z13=p.p3.z-p.p1.z;
	float _r=(_x12*_x13+_y12*_y13+_z12*_z13)/(_x12*_x12+_y12*_y12+_z12*_z12);
	float _xx3=_x13-_r*_x12, _yx3=_y13-_r*_y12, _zx3=_z13-_r*_z12;
	float _x3=sqrt(_xx3*_xx3+_yx3*_yx3+_zx3*_zx3), _12=sqrt(_x12*_x12+_y12*_y12+_z12*_z12);

	float	p2x=0,		p2y=_12,
			p3x=_x3,	p3y=(_x12*_x13+_y12*_y13+_z12*_z13)/_12;
	p.tx2.x=p.txm.a*p2x+p.txm.b*p2y, p.tx2.y=p.txm.c*p2x+p.txm.d*p2y;
	p.tx3.x=p.txm.a*p3x+p.txm.b*p3y, p.tx3.y=p.txm.c*p3x+p.txm.d*p3y;
//	p.tx1.x=		p.tx1.y=0;
//	p.tx2.x=0,	p.tx2.y=_12;
//	p.tx3.x=_x3,	p.tx3.y=(_x12*_x13+_y12*_y13+_z12*_z13)/_12;
}
void			switchMode(Mode *m_)
{
	int ntx_=0, **texture_=0;
	for(;ntx>0;)
	{
		++ntx_;
		texture_=(int**)realloc(texture_, ntx_*sizeof(int*)), texture_[ntx_-1]=m->popTexture();
	}
	m->finish();
	m=m_;
	m->initiate();
	for(int k=ntx_-1;k>=0;k--)
		m->pushTexture(texture_[k]);
	if(ntx_>0)m->enqueueTextures();
}

struct			ContactPoint
{
	bool old;
	Vector3f ra, rb;
	float depth;
	float mass_n, mass_t1, mass_t2, bias;
	float Pn, Pt1, Pt2;
	ContactPoint():old(false), Pn(0), Pt1(0), Pt2(0){}
	ContactPoint(Vector3f const &ra, Vector3f const &rb, float depth):old(false), Pn(0), Pt1(0), Pt2(0), ra(ra), rb(rb), depth(depth){}
//	float lambda_n, lambda_t;
//	ContactPoint():old(false), lambda_n(0), lambda_t(0){}
//	ContactPoint(Vector3f const &ra, Vector3f const &rb, float depth):old(false), lambda_n(0), lambda_t(0), ra(ra), rb(rb), depth(depth){}
};
struct			ContactInfo
{
	int a_idx, b_idx;
	std::vector<ContactPoint> p;
	Vector3f n,//normal pointing out of A
		t1, t2;
	bool persistent;
	ContactInfo():b_idx(-1){}
	ContactInfo(int a_idx, int b_idx, float depth, Vector3f const &ra, Vector3f const &rb, Vector3f const &n, Vector3f const &t1, Vector3f const &t2):a_idx(a_idx), b_idx(b_idx), p(1, ContactPoint(ra, rb, depth)), n(n), t1(t1), t2(t2), persistent(false){}
	void add_cp(Vector3f const &ra, Vector3f const &rb, float depth){p.push_back(ContactPoint(ra, rb, depth));}
};
std::vector<ContactInfo> contacts;
void			initialize()
{
	//int stack_height=2;
	//objects0.resize(stack_height);
	//for(int k=0;k<stack_height;++k)
	//{
	//	auto &A=objects0[k];
	//	A.set_properties_cuboid(150, 200, 50, 0.0002f, 0x453F78);
	//	A.set_position(1000, 50, float(100*(k+1)), 0, 0, 0);
	//	A.set_velocity(0, 0, 0, 0, 0, 0);
	//	A.friction=object_friction;
	//	for(int kt=0, ktEnd=A.tr.size();kt<ktEnd;++kt)
	//		A.tr[kt].tc=rand()<<16|rand();
	//}

	objects0.resize(10);
	for(int k=0, kEnd=objects0.size();k<kEnd;++k)
	{
		auto &A=objects0[k];
		A.set_properties_cuboid(float(rand()%100), float(rand()%100), float(rand()%100), 0.0002f, 0x453F78);//0.0002f
		A.set_position(float(rand()%1000), float(rand()%1000), float(rand()%1000), rand()/1000.f, rand()/1000.f, rand()/1000.f);
		A.set_velocity(0, 0, 0, 0, 0, 0);
		A.friction=object_friction;
		for(int kt=0, ktEnd=A.tr.size();kt<ktEnd;++kt)
			A.tr[kt].tc=rand()<<16|rand();
	}
	if(objects0.size()>1)
	{
		objects0[0].set_properties_cuboid(150, 200, 50, 0.0002f, 0x453F78);//0.2f 0.0002f
		objects0[0].set_position(1000, 200, float(100*(0+1)), 0, 0, 0);
		//objects0[0].set_position(1000, float(50*0), float(100*(0+1)), 0, 0, 0);
			objects0[1].set_properties_cuboid(150, 200, 50, 0.0002f, 0x453F78);
			objects0[1].set_position(1000, 200, float(100*(1+1)), 0, 0, 0);
		//	objects0[1].set_position(1000, float(50*1), float(100*(1+1)), 0, 0, 0);
		for(int kt=0, ktEnd=objects0[0].tr.size();kt<ktEnd;++kt)
		{
			objects0[0].tr[kt].tc=rand()<<16|rand();
			objects0[1].tr[kt].tc=rand()<<16|rand();
		}
	}//*/
/*	for(int k=0, kEnd=objects0.size();k<kEnd;++k)
	{
		auto &A=objects0[k];
		A.set_properties_cuboid(150, 200, 50, 0.0002f, 0x453F78);
	//	A.set_position(1000+float(10*k), float(50*k), float(45*(k+1)), 0, 0, 0);
		A.set_position(1000, float(50*k), float(100*(k+1)), 0, 0, 0);
	//	A.set_position(0, float(50*k), float(100*(k+1)), 0, 0, 0);
		A.set_velocity(0, 0, 0, 0, 0, 0);
	//	A.set_velocity(0, 0, 0, 0.04f, 0.01f, 0.00f);
	//	A.set_velocity(0, 0, 0, 0.01f, 0.02f, 0.03f);
		A.friction=object_friction;
		for(int kt=0, ktEnd=A.tr.size();kt<ktEnd;++kt)
			A.tr[kt].tc=rand()<<16|rand();
	}//*/

	for(int k=0, kEnd=objects0.size();k<kEnd;++k)
	{
		auto &A=objects0[k];
		A.idx=k, A.iterate(0);
	}
	objects=objects0;
}
unsigned		g_message=0;
void			count_active_keys()
{
	kp=kb['W']+kb['A']+kb['S']+kb['D']+kb['T']+kb['G']+kb[VK_UP]+kb[VK_DOWN]+kb[VK_LEFT]+kb[VK_RIGHT]+kb[VK_ADD]+kb[VK_SUBTRACT]+kb[VK_RBUTTON];
}
long			__stdcall WndProc(HWND__ *hWnd, unsigned int message, unsigned int wParam, long lParam)
{
	g_message=message;//
	switch(message)
	{
	case WM_CREATE:
		initialize();
		break;
	case WM_PAINT:
		GetClientRect(hWnd, &R);
		if(h!=R.bottom-R.top||w!=R.right-R.left)
		{
			h=R.bottom-R.top, w=R.right-R.left, centerP.x=X0=w/2, centerP.y=Y0=h/2;
			ClientToScreen(hWnd, &centerP);
			m->resize();
		}
		if(!timer)render();
		break;
	case WM_EXITSIZEMOVE:
		centerP.x=X0, centerP.y=Y0;ClientToScreen(hWnd, &centerP);
		return 0;
	case WM_ACTIVATE:
		if(wParam==WA_INACTIVE)
		{
			memset(kb, 0, 256*sizeof(char));
			kp=0;
			if(timer){KillTimer(hWnd, 0);timer=0;}
		}
		else if(!pause)
			SetTimer(hWnd, 0, TIMER_ELAPSE, 0), timer=1;
		return 0;
	case WM_TIMER:
		if(!pause)
			tick=true;
	//	glParam=lParam, gwParam=wParam;
		render();
		if(pause&&!kp){KillTimer(hWnd, 0);timer=0;}
		return 0;
	case WM_LBUTTONDOWN:
		drag=!drag;
		ShowCursor(!drag);
		if(drag)
		{
			mouseP0.x=short(lParam), mouseP0.y=short(lParam>>16);
			ClientToScreen(hWnd, &mouseP0);
			SetCursorPos(centerP.x, centerP.y);
		}
		else
			SetCursorPos(mouseP0.x, mouseP0.y);
		return 0;
	case WM_MOUSEMOVE://task manager sometimes causes WM_MOUSEMOVE as it updates
		if(drag)
		{
			if(!d_bypass)
			{
				ax+=da_tfov*MOUSE_SENSITIVITY*(X0-short(lParam		)), cax=cos(ax), sax=sin(ax);
				ay+=da_tfov*MOUSE_SENSITIVITY*(Y0-short(lParam>>16	)), cay=cos(ay), say=sin(ay), d_bypass=1;
				for(;ax<0;)ax+=_2pi; for(;ax>_2pi;)ax-=_2pi; for(;ay<0;)ay+=_2pi; for(;ay>_2pi;)ay-=_2pi;
				SetCursorPos(centerP.x, centerP.y);
				if(!timer)render();
			}
			else
				d_bypass=0;
		}
#ifdef AA_MOUSEMOVE_TRIGGER
		else if(!timer)render();//
#endif
		return 0;
	case WM_RBUTTONDOWN:
		kb[VK_RBUTTON]=1, count_active_keys();
		if(!timer){SetTimer(hWnd, 0, TIMER_ELAPSE, 0);timer=1;}
		return 0;
	case WM_RBUTTONUP:
		kb[VK_RBUTTON]=0, count_active_keys();
		return 0;
	case WM_MOUSEWHEEL:
		if(kb[VK_CONTROL])
		{
				 if(short(wParam>>16)>0)	dcam*=2;//wheel up
			else							dcam/=2;//wheel down
		}
		else if(kb[VK_MENU])
		{
				 if(short(wParam>>16)>0)	timescale*=2;
			else							timescale/=2;
		}
		else
		{
				 if(short(wParam>>16)>0)	tanfov/=DTANFOV, fov/=DFOV;
			else							tanfov*=DTANFOV, fov*=DFOV;
			da_tfov=tanfov>1?1:tanfov;
			m->changeFov();
		}
		if(!timer)
			render();
		return 0;
	case WM_KEYDOWN:case WM_SYSKEYDOWN:
		switch(wParam)
		{
		case 'W':case 'A':case 'S':case 'D':case 'T':case 'G':case VK_UP:case VK_DOWN:case VK_LEFT:case VK_RIGHT:case VK_ADD:case VK_SUBTRACT://case 'Q':
			if(!kb[wParam])kb[wParam]=1, count_active_keys();
			if(!timer){SetTimer(hWnd, 0, TIMER_ELAPSE, 0);timer=1;}
			return 0;
		case '6'://previous EPA iteration
			if(n_EPA_iterations)
				EPA_iteration=(EPA_iteration+n_EPA_iterations-1)%n_EPA_iterations;
			if(!timer)
				render();
			break;
		case '7'://next EPA iteration
			if(n_EPA_iterations)
				EPA_iteration=(EPA_iteration+1)%n_EPA_iterations;
			if(!timer)
				render();
			break;
		case VK_RETURN:case 'E'://next iteration
			if(!timer)
			{
				tick=true;
				render();
			}
			break;
		case 'Y'://toggle info
			info=!info;
			if(!timer)
				render();
			break;
		case 'P':case 'F'://toggle pause
			if(timer)
			{
				KillTimer(hWnd, 0);
				timer=false, pause=true;
			}
			else
			{
				timer=true, pause=false;
				SetTimer(hWnd, 0, 10, 0);
			}
			break;
		case VK_F1://planar
			if(mode!=1)
			{
				mode=1, switchMode((Mode*)&lsw);
				kb[wParam]=1;
				if(!timer)render();
			}
			return 0;
		case VK_F2://cylinder panorama
			if(mode!=2&&mode!=7&&mode!=11)
			{
				if(m==(Mode*)&nlsw)
				{
					if(nlsw.rebuffer!=&NonlinearSW::rebuffer1)
					{
						mode=2;
						nlsw.rebuffer=&NonlinearSW::rebuffer1;
						(nlsw.*(nlsw.rebuffer))();
					}
				}
				//else if(m==(Mode*)&nlcl)
				//{
				//	if(nlcl.rebuffer!=0)
				//	{
				//		mode=nlcl.gpu?11:7;
				//		nlcl.rebuffer=0;
				//		nlcl.changeFov();
				//	}
				//}
				else
				{
					mode=2;
					switchMode((Mode*)&nlsw);
					nlsw.rebuffer=&NonlinearSW::rebuffer1;
					(nlsw.*(nlsw.rebuffer))();
				}
				kb[wParam]=1;
				if(!timer)render();
			}
			return 0;
		case VK_F3://vertical cylinder panorama
			if(mode!=3&&mode!=8&&mode!=12)
			{
				if(m==(Mode*)&nlsw)
				{
					if(nlsw.rebuffer!=&NonlinearSW::rebuffer2)
					{
						mode=3;
						nlsw.rebuffer=&NonlinearSW::rebuffer2;
						(nlsw.*(nlsw.rebuffer))();
					}
				}
				//else if(m==(Mode*)&nlcl)
				//{
				//	if(nlcl.rebuffer!=1)
				//	{
				//		mode=nlcl.gpu?12:8;
				//		nlcl.rebuffer=1;
				//		nlcl.changeFov();
				//	}
				//}
				else
				{
					mode=3;
					switchMode((Mode*)&nlsw);
					nlsw.rebuffer=&NonlinearSW::rebuffer2;
					(nlsw.*(nlsw.rebuffer))();
				}
				kb[wParam]=1;
				if(!timer)render();
			}
			return 0;
		case VK_F4://spherical
			if(!kb[VK_MENU]&&mode!=4&&mode!=9&&mode!=13)
			{
				if(m==(Mode*)&nlsw)
				{
					if(nlsw.rebuffer!=&NonlinearSW::rebuffer3)
					{
						mode=4;
						nlsw.rebuffer=&NonlinearSW::rebuffer3;
						(nlsw.*(nlsw.rebuffer))();
					}
				}
				//else if(m==(Mode*)&nlcl)
				//{
				//	if(nlcl.rebuffer!=2)
				//	{
				//		mode=nlcl.gpu?13:9;
				//		nlcl.rebuffer=2;
				//		nlcl.changeFov();
				//	}
				//}
				else
				{
					mode=4;
					switchMode((Mode*)&nlsw);
					nlsw.rebuffer=&NonlinearSW::rebuffer3;
					(nlsw.*(nlsw.rebuffer))();
				}
				kb[wParam]=1;
				if(!timer)render();
			}
			break;
		case VK_F5://fuzzy
			if(mode!=5&&mode!=10&&mode!=14)
			{
				if(m==(Mode*)&nlsw)
				{
					if(nlsw.rebuffer!=&NonlinearSW::rebuffer4)
					{
						mode=5;
						nlsw.rebuffer=&NonlinearSW::rebuffer4;
						(nlsw.*(nlsw.rebuffer))();
					}
				}
				//else if(m==(Mode*)&nlcl)
				//{
				//	if(nlcl.rebuffer!=3)
				//	{
				//		mode=nlcl.gpu?14:10;
				//		nlcl.rebuffer=3;
				//		nlcl.changeFov();
				//	}
				//}
				else
				{
					mode=5;
					switchMode((Mode*)&nlsw);
					nlsw.rebuffer=&NonlinearSW::rebuffer4;
					(nlsw.*(nlsw.rebuffer))();
				}
				kb[wParam]=1;
				if(!timer)render();
			}
			return 0;
		case VK_F6:
			if(mode!=15)
			{
				mode=15, switchMode((Mode*)&psw);
				kb[wParam]=1;
				if(!timer)render();
			}
			return 0;
		case VK_TAB:
			switch(mode)
			{
			case 1:
				mode=6, switchMode((Mode*)&lgl);
				break;
			//case 2:
			//	mode=7, switchMode((Mode*)&nlcl);
			//	nlcl.rebuffer=0;
			//	nlcl.changeFov();
			//	break;
			//case 3:
			//	mode=8, switchMode((Mode*)&nlcl);
			//	nlcl.rebuffer=1;
			//	nlcl.changeFov();
			//	break;
			//case 4:
			//	mode=9, switchMode((Mode*)&nlcl);
			//	nlcl.rebuffer=2;
			//	nlcl.changeFov();
			//	break;
			//case 5:
			//	mode=10, switchMode((Mode*)&nlcl);
			//	nlcl.rebuffer=3;
			//	nlcl.changeFov();
			//	break;

			case 6:
				mode=1, switchMode((Mode*)&lsw);
				break;
			//case 7:
			//	mode=11;
			//	nlcl.switchDevice(1);
			//	break;
			//case 8:
			//	mode=12;
			//	nlcl.switchDevice(1);
			//	break;
			//case 9:
			//	mode=13;
			//	nlcl.switchDevice(1);
			//	break;
			//case 10:
			//	mode=14;
			//	nlcl.switchDevice(1);
			//	break;

			case 11:
				mode=2, switchMode((Mode*)&nlsw);
				nlsw.rebuffer=&NonlinearSW::rebuffer1;
				(nlsw.*(nlsw.rebuffer))();
				break;
			case 12:
				mode=3, switchMode((Mode*)&nlsw);
				nlsw.rebuffer=&NonlinearSW::rebuffer2;
				(nlsw.*(nlsw.rebuffer))();
				break;
			case 13:
				mode=4, switchMode((Mode*)&nlsw);
				nlsw.rebuffer=&NonlinearSW::rebuffer3;
				(nlsw.*(nlsw.rebuffer))();
				break;
			case 14:
				mode=5, switchMode((Mode*)&nlsw);
				nlsw.rebuffer=&NonlinearSW::rebuffer4;
				(nlsw.*(nlsw.rebuffer))();
				break;
			case 15:
				break;
			}
			if(!timer)render();
			return 0;
		case 'R'://reset unpaused
			frame_number=0;
			objects=objects0;
			contacts.resize(0);
			if(pause)//
			{
				timer=true, pause=false;
				SetTimer(hWnd, 0, 10, 0);
			}
			if(!timer)
				tick=true;
			if(!kb[wParam])
			{
				if(kb[VK_CONTROL])
				{
					cam=cam0, dcam=dcam0;
					ax=ax0, ay=ay0, cax=cos(ax), sax=sin(ax), cay=cos(ay), say=sin(ay), da=da0, tanfov=tanfov0, fov=fov0;
					da_tfov=tanfov;
				}
				//if(kb[VK_CONTROL])
				//{
				//	cam.set(0, 0, -1);
				//	ax=0, ay=-.25f*_2pi, cax=1, sax=0, cay=0, say=-1;
				//	kp=0;
				////	memset(kb, 0, 256*sizeof(char));
				//}
				//else
				//{
				//	cam=cam0, dcam=dcam0;
				////	cam.x=camx0, cam.y=camy0, cam.z=camz0, dcam=dcam0;
				//	ax=ax0, ay=ay0, cax=cos(ax), sax=sin(ax), cay=cos(ay), say=sin(ay), da=da0, tanfov=tanfov0, fov=fov0;
				//	da_tfov=tanfov;
				//}
				kb[wParam]=1;
				if(!timer)render();
			}
			return 0;
		case '4'://new scenario
			frame_number=0;
			initialize(), ++scenario_number;
			if(pause)//
			{
				timer=true, pause=false;
				SetTimer(hWnd, 0, 10, 0);
			}
			if(!timer)
			{
				tick=true;
				InvalidateRect(hWnd, 0, false);
			}
			break;
		case 'Q'://reset paused
			frame_number=0;
			objects=objects0;
			contacts.resize(0);
			if(timer)
			{
				KillTimer(hWnd, 0);
				timer=false, pause=true;
			}
		//	tick=true;
			render();
			break;
		case '5'://new scenario paused
			frame_number=0;
			initialize(), ++scenario_number;
			if(timer)
			{
				KillTimer(hWnd, 0);
				timer=false, pause=true;
			}
		//	tick=true;
			render();
			break;
		case 'Z'://turn 180 degrees
			if(!kb[wParam])
			{
				ax+=_2pi/2, ay+=_2pi/2, cax=cos(ax), sax=sin(ax), cay=cos(ay), say=sin(ay);
				kb[wParam]=1;
				if(!timer)render();
			}
			return 0;
		case '0'://random polygons
			{
				m->clearTextures();
				txh=(int*)realloc(txh, sizeof(int)), txh[ntx]=100;
				txw=(int*)realloc(txw, sizeof(int)), txw[ntx]=100;
				int *texture_=(int*)malloc(txh[ntx]*txw[ntx]*sizeof(int));
				for(int k=0;k<txh[ntx]*txw[ntx];++k)texture_[k]=(rand()<<15|rand())&0xFFFFFF;//unsigned char(rand())<<16|unsigned char(rand())<<8|unsigned char(rand());
				m->pushTexture(texture_);
				m->enqueueTextures();
				npols=nRandPols, env=(Triangle0*)realloc(env, npols*sizeof(Triangle0));
				for(int k=0;k<npols;++k)
				{
					auto &p=env[k];
					p.p1.x=float(rand()%500), p.p1.y=float(rand()%500), p.p1.z=float(rand()%500);
					p.p2.x=float(rand()%500), p.p2.y=float(rand()%500), p.p2.z=float(rand()%500);
					p.p3.x=float(rand()%500), p.p3.y=float(rand()%500), p.p3.z=float(rand()%500);
					if(rand()%2)
					{
						float rx=float(rand())/10000, ry=float(rand())/10000, thx=float(rand())/1000, thy=float(rand())/1000;
						p.tx_idx= 0, p.tx1.x=float(rand()%500), p.tx1.y=float(rand()%500), p.tx2.x=float(rand()%500), p.tx2.y=float(rand()%500), p.tx3.x=float(rand()%500), p.tx3.y=float(rand()%500);
						transform(p);
					}
					else
						p.tx_idx=-1, p.color=(rand()<<16|rand())&0xFFFFFF;//unsigned char(rand())<<16|unsigned char(rand())<<8|unsigned char(rand());
				}
				ntpols=nRandTPols, tenv=(Triangle0*)realloc(tenv, ntpols*sizeof(Triangle0));
				for(int k=0;k<ntpols;++k)
				{
					auto &p=tenv[k];
					p.p1.x=float(rand()%500), p.p1.y=float(rand()%500), p.p1.z=float(rand()%500);
					p.p2.x=float(rand()%500), p.p2.y=float(rand()%500), p.p2.z=float(rand()%500);
					p.p3.x=float(rand()%500), p.p3.y=float(rand()%500), p.p3.z=float(rand()%500);
					if(rand()%2)
					{
						float rx=float(rand())/10000, ry=float(rand())/10000, thx=float(rand())/1000, thy=float(rand())/1000;
						p.tx_idx= 0, p.tx1.x=float(rand()%500), p.tx1.y=float(rand()%500), p.tx2.x=float(rand()%500), p.tx2.y=float(rand()%500), p.tx3.x=float(rand()%500), p.tx3.y=float(rand()%500);
						transform(p);
					}
					else
						p.tx_idx=-1, p.color=(rand()<<16|rand())&0xFFFFFF;//unsigned char(rand())<<16|unsigned char(rand())<<8|unsigned char(rand());
				}
			}
			kb[wParam]=1;
			if(!timer)render();
			return 0;
		case '9':
			{
				npols=200, env=(Triangle0*)realloc(env, npols*sizeof(Triangle0));
				int k2, k3, k4=int(sqrt(float(npols)/2));
				for(int k=0;k+1<npols;k+=2)
				{
					auto &p=env[k], &p1=env[k+1];
					k2=(k/2)%k4, k3=(k/2)/k4;
					p .p1.x=1000* k2   +float(rand()%100), p .p1.y=1000* k3   +float(rand()%100), p .p1.z=float(rand()%100);
					p .p2.x=1000*(k2+1)+float(rand()%100), p .p2.y=1000* k3   +float(rand()%100), p .p2.z=float(rand()%100);
					p .p3.x=1000*(k2+1)+float(rand()%100), p .p3.y=1000*(k3+1)+float(rand()%100), p .p3.z=float(rand()%100);
					p1.p1.x=1000*(k2+1)+float(rand()%100), p1.p1.y=1000*(k3+1)+float(rand()%100), p1.p1.z=float(rand()%100);
					p1.p2.x=1000* k2   +float(rand()%100), p1.p2.y=1000*(k3+1)+float(rand()%100), p1.p2.z=float(rand()%100);
					p1.p3.x=1000* k2   +float(rand()%100), p1.p3.y=1000* k3   +float(rand()%100), p1.p3.z=float(rand()%100);
					if(rand()%3&&ntx)
					{
						p.tx_idx=p1.tx_idx=rand()%ntx;
						float rx=float(rand())/1000, ry=float(rand())/1000, thx=float(rand())/1000, thy=float(rand())/1000;
						p .tx1.x=float(rand()%500), p .tx1.y=float(rand()%500), p .txm.a=rx*cos(thx), p .txm.b=ry*cos(thy), p .txm.c=rx*sin(thx), p .txm.d=ry*sin(thy);
						p1.tx1.x=float(rand()%500), p1.tx1.y=float(rand()%500), p1.txm.a=rx*cos(thx), p1.txm.b=ry*cos(thy), p1.txm.c=rx*sin(thx), p1.txm.d=ry*sin(thy);
					}
					else
						p.tx_idx=p1.tx_idx=-1, p.color=p1.color=(rand()<<16|rand())&0xFFFFFF;//unsigned char(rand())<<16|unsigned char(rand())<<8|unsigned char(rand());
				}
			}
			kb[wParam]=1;
			if(!timer)render();
			return 0;
		case '8':
			{
				npols=20, env=(Triangle0*)realloc(env, npols*sizeof(Triangle0));
				int k2, k3, k4=int(sqrt(float(npols)/2));
				for(int k=0;k+1<npols;k+=2)
				{
					auto &p=env[k], &p1=env[k+1];
					k2=(k/2)%k4, k3=(k/2)/k4;
					p .p1.x=float(rand()%1000), p .p1.y=float(rand()%1000), p .p1.z=float(rand()%1000);
					p .p2.x=float(rand()%1000), p .p2.y=float(rand()%1000), p .p2.z=float(rand()%1000);
					p .p3.x=float(rand()%1000), p .p3.y=float(rand()%1000), p .p3.z=float(rand()%1000);
					p1.p1.x=float(rand()%1000), p1.p1.y=float(rand()%1000), p1.p1.z=float(rand()%1000);
					p1.p2.x=float(rand()%1000), p1.p2.y=float(rand()%1000), p1.p2.z=float(rand()%1000);
					p1.p3.x=float(rand()%1000), p1.p3.y=float(rand()%1000), p1.p3.z=float(rand()%1000);
					if(rand()%3&&ntx)
					{
						p.tx_idx=p1.tx_idx=rand()%ntx;
						float rx=float(rand())/1000, ry=float(rand())/1000, thx=float(rand())/1000, thy=float(rand())/1000;
						p .tx1.x=float(rand()%500), p .tx1.y=float(rand()%500), p .txm.a=rx*cos(thx), p .txm.b=ry*cos(thy), p .txm.c=rx*sin(thx), p .txm.d=ry*sin(thy);
						p1.tx1.x=float(rand()%500), p1.tx1.y=float(rand()%500), p1.txm.a=rx*cos(thx), p1.txm.b=ry*cos(thy), p1.txm.c=rx*sin(thx), p1.txm.d=ry*sin(thy);
					}
					else
						p.tx_idx=p1.tx_idx=-1, p.color=p1.color=(rand()<<16|rand())&0xFFFFFF;//unsigned char(rand())<<16|unsigned char(rand())<<8|unsigned char(rand());
				}
			}
			kb[wParam]=1;
			if(!timer)render();
			return 0;
		case '1':
			{
				int k2, k3, k4=int(sqrt(float(npols)/2));if(k4==0)k4=1;
				for(int k=0;k<npols;k+=2)
				{
					auto &p=env[k];
					k2=(k/2)%k4, k3=(k/2)/k4;
									p .p1.x=float(1000* k2   ), p .p1.y=float(1000* k3   ), p .p1.z=float(10*k2*k3);
									p .p2.x=float(1000* k2   ), p .p2.y=float(1000*(k3+1)), p .p2.z=float(10*k2*k3);
									p .p3.x=float(1000*(k2+1)), p .p3.y=float(1000*(k3+1)), p .p3.z=float(10*k2*k3);
									if(p .tx_idx>=0)p .tx1.x=						p .tx1.y=0,						p .tx2.x=(float)txw[p.tx_idx],	p .tx2.y=0,						p .tx3.x=(float)txw[p.tx_idx],	p .tx3.y=(float)txh[p.tx_idx],	transform(p);
					if(k+1<npols){	auto &p1=env[k+1];
									p1.p1.x=float(1000*(k2+1)), p1.p1.y=float(1000*(k3+1)), p1.p1.z=float(10*k2*k3);
									p1.p2.x=float(1000*(k2+1)), p1.p2.y=float(1000* k3   ), p1.p2.z=float(10*k2*k3);
									p1.p3.x=float(1000* k2   ), p1.p3.y=float(1000* k3   ), p1.p3.z=float(10*k2*k3);
									if(p1.tx_idx>=0)p1.tx1.x=(float)txw[p1.tx_idx],	p1.tx1.y=(float)txh[p1.tx_idx],	p1.tx2.x=0,						p1.tx2.y=(float)txh[p1.tx_idx], p1.tx3.x=0,						p1.tx3.y=0,						transform(p1);}
				}
			}
			kb[wParam]=1;
			if(!timer)render();
			return 0;
		case '2':
			for(int k=0;k<npols;++k)
			{
				auto &p=env[k];
				if(p.tx_idx>=0)
				{
					p.tx1.x=p.tx1.y=0;
					p.tx2.x=(float)txw[p.tx_idx], p.tx2.y=0;
					p.tx3.x=(float)txw[p.tx_idx], p.tx3.y=(float)txh[p.tx_idx];
					transform(p);
				}
			}
			kb[wParam]=1;
			if(!timer)render();
			return 0;
		case 'O':
			{
				OpenClipboard(hWnd);
				char *a=(char*)GetClipboardData(CF_OEMTEXT);
				if(!a||strlen(a)<=0){CloseClipboard();return 0;}
				std::string str=a; str+='\\';
				CloseClipboard();
				{
					std::vector<std::wstring> dir;
					std::string str2=str+'*';
					_WIN32_FIND_DATAA data;
					void *hSearch=FindFirstFileA(str2.c_str(), &data);
					for(;;)
					{
						int len=strlen(data.cFileName);
						if(len>=4&&data.cFileName[len-4]=='.'&&(	(data.cFileName[len-3]=='b'||data.cFileName[len-3]=='B')&&(data.cFileName[len-2]=='m'||data.cFileName[len-2]=='M')&&(data.cFileName[len-1]=='p'||data.cFileName[len-1]=='P')
																||	(data.cFileName[len-3]=='p'||data.cFileName[len-3]=='P')&&(data.cFileName[len-2]=='n'||data.cFileName[len-2]=='N')&&(data.cFileName[len-1]=='g'||data.cFileName[len-1]=='G')
																||	(data.cFileName[len-3]=='j'||data.cFileName[len-3]=='J')&&(data.cFileName[len-2]=='p'||data.cFileName[len-2]=='P')&&(data.cFileName[len-1]=='g'||data.cFileName[len-1]=='G')
																))
							dir.push_back(wider(data.cFileName));
						if(!FindNextFileA(hSearch, &data))break;
					}
					FindClose(hSearch);
					std::sort(dir.begin(), dir.end(), [](std::wstring const &a, std::wstring const &b){return StrCmpLogicalW(a.c_str(), b.c_str())==-1;});
					if(dir.size())
					{
						npols=0;
						m->clearTextures();
					}
					Gdiplus::GdiplusStartupInput gdiplusStartupInput;
					ULONG_PTR hgdiplusToken;
					Gdiplus::GdiplusStartup(&hgdiplusToken, &gdiplusStartupInput, 0);
					for(unsigned int k=0;k<dir.size();++k)
					{
						std::wstring str3=wider((char*)str.c_str())+dir[k];
						Gdiplus::Bitmap bitmap(str3.c_str());
						txh=(int*)realloc(txh, (ntx+1)*sizeof(int)), txh[ntx]=bitmap.GetHeight();
						txw=(int*)realloc(txw, (ntx+1)*sizeof(int)), txw[ntx]=bitmap.GetWidth();
						int *texture_=(int*)malloc(txh[ntx]*txw[ntx]*sizeof(int));
						Gdiplus::Rect rt(0, 0, txw[ntx], txh[ntx]);
						Gdiplus::BitmapData data;
						bitmap.LockBits(&rt, Gdiplus::ImageLockModeRead, PixelFormat32bppARGB, &data);
						memcpy(texture_, data.Scan0, txw[ntx]*txh[ntx]*sizeof(int));
						bitmap.UnlockBits(&data);
						m->pushTexture(texture_);
					}
					Gdiplus::GdiplusShutdown(hgdiplusToken);
					if(dir.size())
						m->enqueueTextures();
				}
				str+="!.txt";
				if(loadText(str))
				{
					ntpols=npols=0;
					unsigned char _c0[3]={0};
					for(int i=0, f=0;comtok(str, &i, &f);)
					{
						for(int k=i;k<f;++k)
						{
							switch(str[k])
							{
							case '0':case '1':case '2':case '3':case '4':case '5':case '6':case '7':case '8':case '9':case '`':case '-':
								{
									k=i;
									int prev;
									env=(Triangle0*)realloc(env, ++npols*sizeof(Triangle0));
									memset(env+npols-1, 0, sizeof(Triangle0));
									auto &p=env[npols-1];
									env[npols-1].tx_idx=-1, env[npols-1].color=(rand()<<15|rand())&0xFFFFFF;//unsigned char(rand())<<16|unsigned char(rand())<<8|unsigned char(rand());
									p.p1.x=readDoubleFromBuffer(str, k, f, &k, &prev); if(prev)p.p1.x=npols>1?env[npols-2].p1.x:0;
									if(k!=-1)
									{
										p.p1.y=readDoubleFromBuffer(str, k, f, &k, &prev); if(prev)p.p1.y=npols>1?env[npols-2].p1.y:0;
										if(k!=-1)
										{
											p.p1.z=readDoubleFromBuffer(str, k, f, &k, &prev); if(prev)p.p1.z=npols>1?env[npols-2].p1.z:0;
											if(k!=-1)
											{
												p.p2.x=readDoubleFromBuffer(str, k, f, &k, &prev); if(prev)p.p2.x=npols>1?env[npols-2].p2.x:0;
												if(k!=-1)
												{
													p.p2.y=readDoubleFromBuffer(str, k, f, &k, &prev); if(prev)p.p2.y=npols>1?env[npols-2].p2.y:0;
													if(k!=-1)
													{
														p.p2.z=readDoubleFromBuffer(str, k, f, &k, &prev); if(prev)p.p2.z=npols>1?env[npols-2].p2.z:0;
														if(k!=-1)
														{
															p.p3.x=readDoubleFromBuffer(str, k, f, &k, &prev); if(prev)p.p3.x=npols>1?env[npols-2].p3.x:0;
															if(k!=-1)
															{
																p.p3.y=readDoubleFromBuffer(str, k, f, &k, &prev); if(prev)p.p3.y=npols>1?env[npols-2].p3.y:0;
																if(k!=-1)
																{
																	p.p3.z=readDoubleFromBuffer(str, k, f, &k, &prev); if(prev)p.p3.z=npols>1?env[npols-2].p3.z:0;
																	if(k!=-1)
																	{
																		for(int k2=k;k2<f;++k2)
																		{
																			if(str[k2]=='t'||str[k2]=='T')
																			{
																				if(k2+1<f&&str[k2+1]=='f'||str[k2+1]=='F')
																				{
																					--npols;
																					tenv=(Triangle0*)realloc(tenv, ++ntpols*sizeof(Triangle0));
																					tenv[ntpols-1]=env[npols];

																					tenv[ntpols-1].tx_idx=char(readDoubleFromBuffer(str, k, f, &k, &prev))-1; if(prev)tenv[ntpols-1].tx_idx=ntpols>1?tenv[ntpols-2].tx_idx:0;
																					if(tenv[ntpols-1].tx_idx<0)
																						tenv[ntpols-1].tx_idx=0;
																					if(tenv[ntpols-1].tx_idx>=ntx)
																						tenv[ntpols-1].tx_idx=ntx-1;
																					if(k!=-1)
																					{
																						tenv[ntpols-1].tx1.x=readDoubleFromBuffer(str, k, f, &k, &prev); if(prev)tenv[ntpols-1].tx1.x=ntpols>1?tenv[ntpols-2].tx1.x:0;
																					}
																					if(k!=-1)//... tx1.x\n		read should tell when last value found and when no value found
																					{
																						tenv[ntpols-1].tx1.y=readDoubleFromBuffer(str, k, f, &k, &prev); if(prev)tenv[ntpols-1].tx1.y=ntpols>1?tenv[ntpols-2].tx1.y:0;
																						if(k!=-1)
																						{
																							tenv[ntpols-1].tx2.x=readDoubleFromBuffer(str, k, f, &k, &prev); if(prev)tenv[ntpols-1].tx2.x=ntpols>1?tenv[ntpols-2].tx2.x:0;
																							if(k!=-1)
																							{
																								tenv[ntpols-1].tx2.y=readDoubleFromBuffer(str, k, f, &k, &prev); if(prev)tenv[ntpols-1].tx2.y=ntpols>1?tenv[ntpols-2].tx2.y:0;
																								if(k!=-1)
																								{
																									tenv[ntpols-1].tx3.x=readDoubleFromBuffer(str, k, f, &k, &prev); if(prev)tenv[ntpols-1].tx3.x=ntpols>1?tenv[ntpols-2].tx3.x:0;
																									if(k!=-1)
																									{
																										tenv[ntpols-1].tx3.y=readDoubleFromBuffer(str, k, f, &k, &prev); if(prev)tenv[ntpols-1].tx3.y=ntpols>1?tenv[ntpols-2].tx3.y:0;
																									}
																								}
																							}
																						}
																						transform(tenv[ntpols-1]);
																					}
																					else
																					{
																						tenv[ntpols-1].tx1.x=tenv[ntpols-1].tx1.y=0, tenv[ntpols-1].txm.a=tenv[ntpols-1].txm.d=1, tenv[ntpols-1].txm.b=tenv[ntpols-1].txm.c=0;
																						antitransform_transparent(ntpols-1);
																					}
																				}
																				else
																				{
																					env[npols-1].tx_idx=char(readDoubleFromBuffer(str, k, f, &k, &prev))-1; if(prev)env[npols-1].tx_idx=npols>1?env[npols-2].tx_idx:0;
																					if(env[npols-1].tx_idx<0)
																						env[npols-1].tx_idx=0;
																					if(env[npols-1].tx_idx>=ntx)
																						env[npols-1].tx_idx=ntx-1;
																					if(k!=-1)
																					{
																						env[npols-1].tx1.x=readDoubleFromBuffer(str, k, f, &k, &prev); if(prev)env[npols-1].tx1.x=npols>1?env[npols-2].tx1.x:0;
																					}
																					if(k!=-1)//... tx1.x\n		read should tell when last value found and when no value found
																					{
																						env[npols-1].tx1.y=readDoubleFromBuffer(str, k, f, &k, &prev); if(prev)env[npols-1].tx1.y=npols>1?env[npols-2].tx1.y:0;
																						if(k!=-1)
																						{
																							env[npols-1].tx2.x=readDoubleFromBuffer(str, k, f, &k, &prev); if(prev)env[npols-1].tx2.x=npols>1?env[npols-2].tx2.x:0;
																							if(k!=-1)
																							{
																								env[npols-1].tx2.y=readDoubleFromBuffer(str, k, f, &k, &prev); if(prev)env[npols-1].tx2.y=npols>1?env[npols-2].tx2.y:0;
																								if(k!=-1)
																								{
																									env[npols-1].tx3.x=readDoubleFromBuffer(str, k, f, &k, &prev); if(prev)env[npols-1].tx3.x=npols>1?env[npols-2].tx3.x:0;
																									if(k!=-1)
																									{
																										env[npols-1].tx3.y=readDoubleFromBuffer(str, k, f, &k, &prev); if(prev)env[npols-1].tx3.y=npols>1?env[npols-2].tx3.y:0;
																									}
																								}
																							}
																						}
																						transform(env[npols-1]);
																					}
																					else
																					{
																						env[npols-1].tx1.x=env[npols-1].tx1.y=0, env[npols-1].txm.a=env[npols-1].txm.d=1, env[npols-1].txm.b=env[npols-1].txm.c=0;
																						antitransform(npols-1);
																					}
																				}
																				break;
																			}
																			else if(str[k2]=='c'||str[k2]=='C')
																			{
																				if(k2+1<f&&str[k2+1]=='f'||str[k2+1]=='F')
																				{
																					auto _c=(unsigned char*)&env[npols-1].color;
																				//	unsigned char _c[3]={0};
																					if(k!=-1)
																					{
																						_c[2]=unsigned char(readDoubleFromBuffer(str, k2+2, f, &k, &prev)); if(prev)_c[2]=npols>1?_c0[2]:0;//r
																						if(k!=-1)
																						{
																							_c[1]=unsigned char(readDoubleFromBuffer(str, k, f, &k, &prev)); if(prev)_c[1]=npols>1?_c0[1]:0;//g
																							if(k!=-1)
																							{
																								_c[0]=unsigned char(readDoubleFromBuffer(str, k, f, &k, &prev)); if(prev)_c[0]=npols>1?_c0[0]:0;//b
																							}
																						}
																					}
																					else
																						_c[0]=rand(), _c[1]=rand(), _c[2]=rand();
																					memcpy(_c0, _c, 3*sizeof(char));
																				}
																				else
																				{
																					auto _c=(unsigned char*)&env[npols-1].color;
																				//	unsigned char _c[3]={0};
																					if(k!=-1)
																					{
																						_c[2]=unsigned char(readDoubleFromBuffer(str, k2+1, f, &k, &prev)); if(prev)_c[2]=npols>1?_c0[2]:0;//r
																						if(k!=-1)
																						{
																							_c[1]=unsigned char(readDoubleFromBuffer(str, k, f, &k, &prev)); if(prev)_c[1]=npols>1?_c0[1]:0;//g
																							if(k!=-1)
																							{
																								_c[0]=unsigned char(readDoubleFromBuffer(str, k, f, &k, &prev)); if(prev)_c[0]=npols>1?_c0[0]:0;//b
																							}
																						}
																					}
																					else
																						_c[0]=rand(), _c[1]=rand(), _c[2]=rand();
																				//	env[npols-1].color=_c[0]<<16|_c[1]<<8|_c[2];
																					memcpy(_c0, _c, 3*sizeof(char));
																				}
																				break;
																			}
																		}
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
								break;
							default:
								continue;
							}
							break;
						}
					}
				}
			}
			kb[wParam]=1;
			if(!timer)render();
			return 0;
		case 'X':
			PostQuitMessage(0);
			return 0;
		}
		kb[wParam]=1;
		break;
	case WM_KEYUP:
		kb[wParam]=0;
		switch(wParam)
		{
		case 'W':case 'A':case 'S':case 'D':case 'T':case 'G':case VK_UP:case VK_DOWN:case VK_LEFT:case VK_RIGHT:case VK_ADD:case VK_SUBTRACT://case 'Q':
			count_active_keys();
			//if(kp>0)kp--;
			break;
		}
		return 0;
	case WM_CLOSE:
		PostQuitMessage(0);
		return 0;
	}
	return DefWindowProcA(hWnd, message, wParam, lParam);
}
int				__stdcall WinMain(HINSTANCE__ *hInstance, HINSTANCE__*, char*, int nCmdShow)
{
	tagWNDCLASSEXA wndClassEx={sizeof(tagWNDCLASSEXA), CS_HREDRAW|CS_VREDRAW|CS_DBLCLKS, WndProc, 0, 0, hInstance, LoadIconA(0, (char*)0x00007F00), LoadCursorA(0, (char*)0x00007F00), (HBRUSH__*)(COLOR_WINDOW+1), 0, "New format", 0};
	RegisterClassExA(&wndClassEx);
	ghWnd=CreateWindowExA(0, wndClassEx.lpszClassName, "Physics 3/AWM20190805", WS_CAPTION|WS_SYSMENU|WS_THICKFRAME|WS_MINIMIZEBOX|WS_MAXIMIZEBOX|WS_CLIPCHILDREN, CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, 0, 0, hInstance, 0);
	ShowWindow(ghWnd, nCmdShow);
	
		QueryPerformanceCounter(&li), srand(li.LowPart);
		nticks=li.QuadPart;
		m=(Mode*)&lgl;
	//	m=(Mode*)&lsw;
	//	nlsw.rebuffer=&NonlinearSW::rebuffer3;
	//	nlcl.rebuffer=2;
		ghDC=GetDC(ghWnd), m->initiate();

	tagMSG msg;
	for(;GetMessageA(&msg, 0, 0, 0);)TranslateMessage(&msg), DispatchMessageA(&msg);
	
		m->finish(), ReleaseDC(ghWnd, ghDC);
		free(env), free(txh), free(txw);

	return msg.wParam;
}

void			draw_projected_xz(PhysObject const &A, int x0, int y0, int color)
{
	for(int kt=0, ktEnd=A.tr.size();kt<ktEnd;++kt)
	{
		auto &face=A.tr[kt];
		auto &v1=A.w[face.a], &v2=A.w[face.b], &v3=A.w[face.c];
		m->draw_line(x0+v1.x, y0-v1.z, x0+v2.x, y0-v2.z, color);
		m->draw_line(x0+v2.x, y0-v2.z, x0+v3.x, y0-v3.z, color);
		m->draw_line(x0+v3.x, y0-v3.z, x0+v1.x, y0-v1.z, color);
	}
}
void			draw_line_projected_xz(Vector3f const &p1, Vector3f const &p2, int x0, int y0, int lineColor){m->draw_line(x0+p1.x, y0-p1.z, x0+p2.x, y0-p2.z, lineColor);}
void			draw_projected_yz(PhysObject const &A, int x0, int y0, int color)
{
	for(int kt=0, ktEnd=A.tr.size();kt<ktEnd;++kt)
	{
		auto &face=A.tr[kt];
		auto &v1=A.w[face.a], &v2=A.w[face.b], &v3=A.w[face.c];
		m->draw_line(x0+v1.y, y0-v1.z, x0+v2.y, y0-v2.z, color);
		m->draw_line(x0+v2.y, y0-v2.z, x0+v3.y, y0-v3.z, color);
		m->draw_line(x0+v3.y, y0-v3.z, x0+v1.y, y0-v1.z, color);
	}
}
void			draw_line_projected_yz(Vector3f const &p1, Vector3f const &p2, int x0, int y0, int lineColor){m->draw_line(x0+p1.y, y0-p1.z, x0+p2.y, y0-p2.z, lineColor);}
void			draw_projected_xy(PhysObject const &A, int x0, int y0, int color)
{
	for(int kt=0, ktEnd=A.tr.size();kt<ktEnd;++kt)
	{
		auto &face=A.tr[kt];
		auto &v1=A.w[face.a], &v2=A.w[face.b], &v3=A.w[face.c];
		m->draw_line(x0+v1.x, y0-v1.y, x0+v2.x, y0-v2.y, color);
		m->draw_line(x0+v2.x, y0-v2.y, x0+v3.x, y0-v3.y, color);
		m->draw_line(x0+v3.x, y0-v3.y, x0+v1.x, y0-v1.y, color);
	}
}
void			draw_line_projected_xy(Vector3f const &p1, Vector3f const &p2, int x0, int y0, int lineColor){m->draw_line(x0+p1.x, y0-p1.y, x0+p2.x, y0-p2.y, lineColor);}
void			draw_3projections(PhysObject const &A, int color)
{
	int x1_4=X0/2, y1_4=Y0/2, x3_4=X0+x1_4, y3_4=Y0+y1_4;
	draw_projected_xz(A, x1_4, Y0, color);
	draw_projected_yz(A, x3_4, Y0, color);
//	draw_projected_xy(A, x1_4, y3_4, color);
	draw_projected_xy(A, x1_4, h, color);
}
void			draw_line_3projections(Vector3f const &p1, Vector3f const &p2, int color)
{
	int x1_4=X0/2, y1_4=Y0/2, x3_4=X0+x1_4, y3_4=Y0+y1_4;
	draw_line_projected_xz(p1, p2, x1_4, Y0, color);
	draw_line_projected_yz(p1, p2, x3_4, Y0, color);
//	draw_line_projected_xy(p1, p2, x1_4, y3_4, color);
	draw_line_projected_xy(p1, p2, x1_4, h, color);
}
//void			draw_projected_yz(PhysObject const &A, int color)
//{
//	for(int kt=0, ktEnd=A.tr.size();kt<ktEnd;++kt)
//	{
//		auto &face=A.tr[kt];
//		auto &v1=A.w[face.a], &v2=A.w[face.b], &v3=A.w[face.c];
//		m->draw_line(X0+v1.y, Y0-v1.z, X0+v2.y, Y0-v2.z, color);
//		m->draw_line(X0+v2.y, Y0-v2.z, X0+v3.y, Y0-v3.z, color);
//		m->draw_line(X0+v3.y, Y0-v3.z, X0+v1.y, Y0-v1.z, color);
//	}
//	//for(int kv=0, kvEnd=A.w.size();kv<kvEnd;++kv)
//	//{
//	//	auto &v0=A.w[kv], &v1=A.w[(kv+1)%A.w.size()];
//	//	float x1=X0+v0.y, y1=Y0-v0.z, x2=X0+v1.y, y2=Y0-v1.z;
//	////	float x1=X0+v0.y*0.5f, y1=Y0+v0.z*0.5f, x2=X0+v1.y*0.5f, y2=Y0+v1.z*0.5f;
//	//	m->draw_line(x1, y1, x2, y2, 0);
//	//}
//}
//void			draw_line_projected_yz(Vector3f const &p1, Vector3f const &p2, int lineColor){m->draw_line(X0+p1.y, Y0-p1.z, X0+p2.y, Y0-p2.z, lineColor);}

void			draw_minkowski_difference(PhysObject const &A, PhysObject const &B, int color)
{
	//Vector3f zero(0, 0, 0);
	//std::vector<qh_vertex> m_diff(A.w.size()*B.w.size());
	//auto it=m_diff.begin();
	for(int ka=0, kaEnd=A.w.size();ka<kaEnd;++ka)
	{
		auto &pa=A.w[ka];
		//for(int kb=0, kbEnd=B.w.size();kb<kbEnd;++kb)
		//{
		//	auto &pb=B.w[kb];
		//	Vector3f p=pa-pb;
		//	it->x=p.x, it->y=p.y, it->z=p.z;
		//	++it;
		//}
		Vector3f p1=pa-B.w[0];
		for(int kb=1, kbEnd=B.w.size();kb<kbEnd;++kb)
		{
			auto &pb=B.w[kb];
			Vector3f p2=pa-pb;
			m->draw_line(p1, p2, color);
			p1=p2;
		}
		//for(int kb=0, kbEnd=B.w.size();kb<kbEnd;++kb)
		//{
		//	auto &pb=B.w[kb];
		//	Vector3f p=pb-pa;
		//	m->draw_point(p, color);
		//	p.x+=1;
		//	m->draw_point(p, color);
		//	p.y+=1;
		//	m->draw_point(p, color);
		//	p.x-=1;
		//	m->draw_point(p, color);
		//	//m->draw_line(zero, p, color);
		//}
	}
	//qh_mesh hull=qh_quickhull3d(&m_diff[0], m_diff.size());
	//auto &v1=hull.vertices[hull.indices[0]];
	//Vector3f p1(v1.x, v1.y, v1.z);
	//for(unsigned k=1;k<hull.nindices;++k)
	//{
	//	auto &v2=hull.vertices[hull.indices[k]];
	//	Vector3f p2(v2.x, v2.y, v2.z);
	//	m->draw_line(p1, p2, 0xFF0000);
	//	p1=p2;
	//}
	//qh_free_mesh(hull);
}
Vector3f		g_ca, g_cb;
Vector3f		vector_to_origin;

bool			inclusion_point_convex_polyhedron(PhysObject const &A, Vector3f const &p)
{
	for(int kt=0, ktEnd=A.tr.size();kt<ktEnd;++kt)
	{
		auto &tr=A.tr[kt];
		const Vector3f &v1=A.w[tr.a], &v2=A.w[tr.b], &v3=A.w[tr.c];
		Vector3f normal=(v2-v1).cross(v3-v1);//normal pointing out of A
		if(normal.dot(p-v1)>0)
		{
			//for(int ke=0, keEnd=A.e.size();ke<keEnd;++ke)
			//{
			//	auto &e=A.e[ke];
			//	m->draw_line(A.w[e.v1], A.w[e.v2], 0);
			//}
			//m->print(p, "outside");
			return false;
		}
	}
	return true;
}

struct			PolytopeVertex
{
	Vector3f p;
	int a_idx, b_idx;
	PolytopeVertex():a_idx(-1), b_idx(-1){}
	PolytopeVertex(Vector3f const &p, int a_idx, int b_idx):p(p), a_idx(a_idx), b_idx(b_idx){}
	void set(Vector3f const &p, int a_idx, int b_idx){this->p=p, this->a_idx=a_idx, this->b_idx=b_idx;}
};
struct			Simplex
{
	PolytopeVertex a, b, c, d;
	int count;
	Simplex():count(0){}
	Simplex(PolytopeVertex const &a):count(1), a(a){}
	Simplex(PolytopeVertex const &a, PolytopeVertex const &b):count(2), a(a), b(b){}
	Simplex(PolytopeVertex const &a, PolytopeVertex const &b, PolytopeVertex const &c):count(3), a(a), b(b), c(c){}
	Simplex(Vector3f const &a, int a_idx, int b_idx):count(1), a(a, a_idx, b_idx){}
	void add(Vector3f const &p, int a_idx, int b_idx)
	{
		if(count<4)
			(&a)[count].set(p, a_idx, b_idx), ++count;
	}
};
struct			Face
{
	PolytopeVertex v1, v2, v3;
	Vector3f n;
	float origin_distance;
	bool contains_origin;
	Face():origin_distance(infinity){}
	Face(PolytopeVertex const &v1, PolytopeVertex const &v2, PolytopeVertex const &v3):v1(v1), v2(v2), v3(v3)
	{
		n=(v2.p-v1.p).cross(v3.p-v1.p);
		n/=n.magnitude();
		origin_distance=v1.p.dot(n);
		if(origin_distance<0)
		{
			PolytopeVertex temp=this->v2; this->v2=this->v3, this->v3=temp;
			origin_distance=-origin_distance, n=-n;
		}
	}
	void set(PolytopeVertex const &v1, PolytopeVertex const &v2, PolytopeVertex const &v3, bool contains_origin)
	{
		this->v1=v1, this->v2=v2, this->v3=v3, this->contains_origin=contains_origin;
		n=(v2.p-v1.p).cross(v3.p-v1.p);
		float mag=n.magnitude();
		if(mag)
			n/=mag;
		//else
		//	int LOL_1=0;
		origin_distance=v1.p.dot(n);
		if(origin_distance<0)
		{
			PolytopeVertex temp=this->v2; this->v2=this->v3, this->v3=temp;
			origin_distance=-origin_distance, n=-n;
		}
	}
};
std::vector<std::vector<Face>> g_faces;
int g_a_idx=0,
	g_b_idx=6;//1
int				g_c=-1;
Vector3f		support(PhysObject const &S, Vector3f const &d, int &index)
{
	index=0;
	float max_dot=S.w[0].dot(d);
	for(int k=1, kEnd=S.w.size();k<kEnd;++k)
	{
		float dot=S.w[k].dot(d);
		if(max_dot<dot)
			max_dot=dot, index=k;
	}
	return S.w[index];
}
void			get_uv_line(Vector3f const &a, Vector3f const &b, float &u, float &v, Vector3f &t)
{
	t=b-a;
	t/=t.mag_sq();
	v=-a.dot(t), u=b.dot(t);
}
void			get_uv_line(Vector3f const &a, Vector3f const &b, float &u, float &v)
{
	Vector3f t=b-a;
	t/=t.mag_sq();
	v=-a.dot(t), u=b.dot(t);
}
void			get_uvw_face_sign_only(Vector3f const &a, Vector3f const &b, Vector3f const &c, float &u, float &v, float &w)
{
	Vector3f n=(b-a).cross(c-a);
	u=b.cross(c).dot(n);
	v=c.cross(a).dot(n);
	w=a.cross(b).dot(n);
}
void			get_uvw_face(Vector3f const &a, Vector3f const &b, Vector3f const &c, float &u, float &v, float &w)
{//https://math.stackexchange.com/questions/128991/how-to-calculate-the-area-of-a-3d-triangle
	Vector3f n=(b-a).cross(c-a);
	float inv_area_sq=1/n.mag_sq();//double area of abc squared
	u=b.cross(c).dot(n)*inv_area_sq;
	v=c.cross(a).dot(n)*inv_area_sq;
	w=a.cross(b).dot(n)*inv_area_sq;
}
inline float	tetrahedron_volume(Vector3f const &ab, Vector3f const &ac, Vector3f const &ad){return ab.cross(ac).dot(ad);}
void			set_GJK_direction_line(Vector3f const &t, Vector3f const &b, Vector3f &d)//set d orthogonal to t, pointing in general direction of the origin from p
{
	d=t.triple_product(-b, t);
}
void			set_GJK_direction_face(Vector3f const &tAB, Vector3f const &tAC, Vector3f const &a, Vector3f &d)
{
	d=tAB.cross(tAC);
	if(d.dot(a)>0)
		d=-d;
}
bool			intersection_polyhedra_GJK(PhysObject const &P, PhysObject const &Q, Vector3f const &initial_axis, Simplex &s)//initial axis: any starting direction works (relative velocity?), returns closest simplex to origin
{//https://gist.github.com/vurtun/29727217c269a2fbf4c0ed9a1d11cb40
	int p_idx, q_idx;
	Vector3f a=support(P, initial_axis, p_idx)-support(Q, -initial_axis, q_idx);
	s=Simplex(a, p_idx, q_idx);
	//Simplex s(a);
	Vector3f d=-a;
	bool contains_origin=false;
	for(;;)
	{
		a=support(P, d, p_idx)-support(Q, -d, q_idx);
		for(int k=0;k<s.count;++k)//terminate if spupport returns duplicate results
			if(p_idx==(&s.a)[k].a_idx&&q_idx==(&s.a)[k].b_idx)
				return false;
		if(a.dot(d)<0)
			return false;
		s.add(a, p_idx, q_idx);
		//nearest simplex (s) -> s,d,contains_origin
		switch(s.count)
		{
		case 1:
			contains_origin=!s.a.p.mag_sq();
			break;
		case 2://line segment: find nearest point or line to the origin
			{
				Vector3f t;
				float u, v;
				get_uv_line(s.a.p, s.b.p, u, v, t);
				Vector3f p;
				if(u<=0)
					s=Simplex(s.b), d=-s.b.p, p=s.b.p;
				else if(v<=0)//should never hit (a is old)
					s=Simplex(s.a), d=-s.a.p, p=s.a.p;
				else
				{
					p=u*s.a.p+v*s.b.p;
					set_GJK_direction_line(t, s.b.p, d);
				}
				contains_origin=!p.mag_sq();
			}
			break;
		case 3://triangle
			{
				if(P.idx==g_a_idx&&Q.idx==g_b_idx)//
					g_faces.push_back(std::vector<Face>(1, Face(s.a, s.b, s.c))), ++n_GJK_iterations, ++n_EPA_iterations;//
				//Vector3f tBC, tCA;
				//float uBC, vBC, uCA, vCA;
				Vector3f p;
				Vector3f tAB, tBC, tCA;
				float uAB, vAB, uBC, vBC, uCA, vCA;
				get_uv_line(s.a.p, s.b.p, uAB, vAB, tAB);//ab is old
				get_uv_line(s.b.p, s.c.p, uBC, vBC, tBC);
				get_uv_line(s.c.p, s.a.p, uCA, vCA, tCA);
				if(vAB<=0&&uCA<=0)//a is old
					s=Simplex(s.a), d=-s.a.p, p=s.a.p;
				else if(uAB<=0&&vBC<=0)//b is old
					s=Simplex(s.b), d=-s.b.p, p=s.b.p;
				else if(uBC<=0&&vCA<=0)
					s=Simplex(s.c), d=-s.c.p, p=s.c.p;
				else
				{
					float uABC, vABC, wABC;
					get_uvw_face(s.a.p, s.b.p, s.c.p, uABC, vABC, wABC);
					//float area=(s.b.p-s.a.p).cross(s.c.p-s.a.p).magnitude(),//ABC
					//	uABC=s.b.p.cross(s.c.p).magnitude()/area,//OBC/ABC
					//	vABC=s.c.p.cross(s.a.p).magnitude()/area;//OCA/ABC
					////	wABC=s.a.p.cross(s.b.p).magnitude()/area;//OAB/ABC
					if(uAB>0&&vAB>0&&wABC<=0)//should never hit, ab is old
						s=Simplex(s.a, s.b), set_GJK_direction_line(tAB, s.b.p, d), p=uAB*s.a.p+vAB*s.b.p;
					else if(uBC>0&&vBC>0&&uABC<=0)
						s=Simplex(s.b, s.c), set_GJK_direction_line(tBC, s.c.p, d), p=uBC*s.b.p+vBC*s.c.p;
					else if(uCA>0&&vCA>0&&vABC<=0)
						s=Simplex(s.a, s.c), set_GJK_direction_line(tCA, s.a.p, d), p=uCA*s.b.p+vCA*s.c.p;
					else
					{
						set_GJK_direction_face(s.b.p-s.a.p, -tCA, s.a.p, d);//the origin is inside ABC prism
						p=uABC*s.a.p+vABC*s.b.p+wABC*s.c.p;
					}
				}
				contains_origin=p.mag_sq()<1e-6;
			}
			break;
		case 4://tetrahedron
			{
				if(P.idx==g_a_idx&&Q.idx==g_b_idx)//
				{
					g_faces.push_back(std::vector<Face>()), ++n_GJK_iterations, ++n_EPA_iterations;
					g_faces.rbegin()->push_back(Face(s.a, s.b, s.c));//
					g_faces.rbegin()->push_back(Face(s.a, s.b, s.d));//
					g_faces.rbegin()->push_back(Face(s.b, s.c, s.d));//
					g_faces.rbegin()->push_back(Face(s.a, s.c, s.d));//
				}
				Vector3f p;
				Vector3f tAB, tBC, tCA, tBD, tDC, tAD;//6 edges
				float uAB, vAB,  uBC, vBC,  uCA, vCA,  uBD, vBD,  uDC, vDC,  uAD, vAD;
				get_uv_line(s.a.p, s.b.p, uAB, vAB, tAB);
				get_uv_line(s.b.p, s.c.p, uBC, vBC, tBC);
				get_uv_line(s.c.p, s.a.p, uCA, vCA, tCA);
				get_uv_line(s.b.p, s.d.p, uBD, vBD, tBD);
				get_uv_line(s.d.p, s.c.p, uDC, vDC, tDC);
				get_uv_line(s.a.p, s.d.p, uAD, vAD, tAD);
				if(vAB<=0&&uCA<=0&&vAD<=0)//a is old
					s=Simplex(s.a), d=-s.a.p, p=s.a.p;
				else if(uAB<=0&&vBC<=0&&vBD<=0)//b is old
					s=Simplex(s.b), d=-s.b.p, p=s.b.p;
				else if(uBC<=0&&vCA<=0&&uDC<=0)//c is old
					s=Simplex(s.c), d=-s.c.p, p=s.c.p;
				else if(uBD<=0&&vDC<=0&&uAD<=0)
					s=Simplex(s.d), d=-s.d.p, p=s.d.p;
				else
				{
					float uADB, vADB, wADB,  uACD, vACD, wACD,  uCBD, vCBD, wCBD,  uABC, vABC, wABC;
					get_uvw_face(s.a.p, s.d.p, s.b.p, uADB, vADB, wADB);
					get_uvw_face(s.a.p, s.c.p, s.d.p, uACD, vACD, wACD);
					get_uvw_face(s.c.p, s.b.p, s.d.p, uCBD, vCBD, wCBD);
					get_uvw_face(s.a.p, s.b.p, s.c.p, uABC, vABC, wABC);
					if(wABC<=0&&vADB<=0&&uAB>0&&vAB>0)//ab is old
						s=Simplex(s.a, s.b), set_GJK_direction_line(tAB, s.b.p, d), p=uAB*s.a.p+vAB*s.b.p;
					else if(uABC<=0&&wCBD<=0&&uBC>0&&vBC>0)//bc is old
						s=Simplex(s.b, s.c), set_GJK_direction_line(tBC, s.c.p, d), p=uBC*s.b.p+vBC*s.c.p;
					else if(vABC<=0&&wACD<=0&&uCA>0&&vCA>0)//ca is old
						s=Simplex(s.a, s.c), set_GJK_direction_line(tCA, s.a.p, d), p=uCA*s.c.p+vCA*s.a.p;
					else if(vCBD<=0&&uACD<=0&&uDC>0&&vDC>0)
						s=Simplex(s.c, s.d), set_GJK_direction_line(tDC, s.c.p, d), p=uDC*s.d.p+vDC*s.c.p;
					else if(vACD<=0&&wADB<=0&&uAD>0&&vAD>0)
						s=Simplex(s.a, s.d), set_GJK_direction_line(tAD, s.d.p, d), p=uAD*s.a.p+vAD*s.d.p;
					else if(uCBD<=0&&uADB<=0&&uBD>0&&vBD>0)
						s=Simplex(s.b, s.d), set_GJK_direction_line(tBD, s.d.p, d), p=uBD*s.b.p+vBD*s.d.p;
					else
					{
						float inv_vol=tetrahedron_volume(s.b.p-s.a.p, s.c.p-s.a.p, s.d.p-s.a.p);
						inv_vol=inv_vol?1/inv_vol:1;
						float uABCD=tetrahedron_volume(s.c.p, s.d.p, s.b.p)*inv_vol;
						float vABCD=tetrahedron_volume(s.c.p, s.a.p, s.d.p)*inv_vol;
						float wABCD=tetrahedron_volume(s.d.p, s.a.p, s.b.p)*inv_vol;
						float xABCD=tetrahedron_volume(s.b.p, s.a.p, s.c.p)*inv_vol;
						if(xABCD<=0&&uABC>0&&vABC>0&&wABC>0)//abc is old
							s=Simplex(s.a, s.b, s.c), set_GJK_direction_face(tAB, -tCA, s.a.p, d), p=uABC*s.a.p+vABC*s.b.p+wABC*s.c.p;
						else if(uABCD<=0&&uCBD>0&&vCBD>0&&wCBD>0)
							s=Simplex(s.c, s.b, s.d), set_GJK_direction_face(-tBC, -tDC, s.c.p, d), p=uCBD*s.c.p+vCBD*s.b.p+wCBD*s.d.p;
						else if(vABCD<=0&&uACD>0&&vACD>0&&wACD>0)
							s=Simplex(s.a, s.c, s.d), set_GJK_direction_face(-tCA, tAD, s.a.p, d), p=uACD*s.a.p+vACD*s.c.p+wACD*s.d.p;
						else if(wABCD<=0&&uADB>0&&vADB>0&&wADB>0)
							s=Simplex(s.a, s.b, s.d), set_GJK_direction_face(tAB, tAD, s.a.p, d), p=uADB*s.a.p+vADB*s.d.p+wADB*s.b.p;
						else
							return true;//your zipper is open
					}
				}
				contains_origin=p.mag_sq()<1e-6;
			}
			break;
		}
		if(contains_origin)
			break;
	}
	return contains_origin;
}
void			make_orthogonal(Vector3f &v)
{
	float xy=sqrt(v.x*v.x+v.y*v.y);
	if(xy)
		v.set(-v.z*v.x/xy, -v.z*v.y/xy, xy);
	else//v is vertical
		v.set(1, 0, 0);
}
//Vector3f		orthogonal_vector(Vector3f const &v)
//{
//	float xy=sqrt(v.x*v.x+v.y*v.y);
//	if(xy)
//		return Vector3f(-v.z*v.x/xy, -v.z*v.y/xy, xy);
//	return Vector3f(1, 0, 0);//v is vertical
//}
bool			intersection_EPA(PhysObject const &A, PhysObject const &B, Vector3f &normal, float &depth, Vector3f &ca, Vector3f &cb)
{	//https://hero.handmade.network/forums/code-discussion/t/1400-expanding_polytope_algorithm_in_3d_confusion.
	//https://github.com/CG-F15-9-Rutgers/SteerLite/blob/master/steerlib/src/GJK_EPA.cpp
	//https://github.com/kevinmoran/GJK/blob/master/GJK.h
	Simplex s;
	if(A.idx==g_a_idx&&B.idx==g_b_idx)//
		g_faces.resize(0), n_GJK_iterations=0, n_EPA_iterations=0;//
	//if(frame_number==32)
	//if(frame_number==35)
	//if(frame_number==32)
	if(frame_number==574)
		int LOL_1=0;
	if(intersection_polyhedra_GJK(A, B, A.p-B.p, s))
	{
		int a_idx, b_idx, a_idx2, b_idx2;
		float u, v, w;
		int s0_count=s.count;
		std::vector<Face> faces;
		switch(s.count)//expand the simplex to tetrahedron
		{
		case 1://point
			{
				Vector3f point=support(A, s.a.p, a_idx)-support(B, -s.a.p, b_idx);
				if(a_idx==s.a.a_idx&&b_idx==s.a.b_idx)
				{
					depth=point.magnitude();
					normal=-point;
					normal*=1/depth;
					ca=A.w[s.a.a_idx];
					cb=B.w[s.a.b_idx];
					return true;
				}
				s.add(point, a_idx, b_idx);
			}//continued in line case
		case 2://line segment contains origin
			{
				Vector3f edge=s.b.p-s.a.p, n=edge.triple_product(s.b.p, edge);
				Vector3f point=support(A, n, a_idx)-support(B, -n, b_idx);
				Vector3f point2=support(A, -n, a_idx2)-support(B, n, b_idx2);
				if(a_idx==s.a.a_idx&&b_idx==s.a.b_idx||a_idx==s.b.a_idx&&b_idx==s.b.b_idx)//farthest point is already in simplex
				{
					depth=point.magnitude();
					//normal=s.b.p-s.a.p;
					//normal*=1/normal.magnitude();
					//make_orthogonal(normal);
					normal=-point;
					normal*=1/depth;

					get_uv_line(s.a.p, s.b.p, u, v);
					ca=u*A.w[s.a.a_idx]+v*A.w[s.b.a_idx];
					cb=u*B.w[s.a.b_idx]+v*B.w[s.b.b_idx];
					return true;
				}
				else if(a_idx2==s.a.a_idx&&b_idx2==s.a.b_idx||a_idx2==s.b.a_idx&&b_idx2==s.b.b_idx)
				{
					depth=point2.magnitude();
					//normal=s.b.p-s.a.p;
					//normal*=1/normal.magnitude();
					//make_orthogonal(normal);
					normal=-point2;
					normal*=1/depth;

					get_uv_line(s.a.p, s.b.p, u, v);
					ca=u*A.w[s.a.a_idx]+v*A.w[s.b.a_idx];
					cb=u*B.w[s.a.b_idx]+v*B.w[s.b.b_idx];
					return true;
				}
				s.add(point, a_idx, b_idx);
				s.add(point2, a_idx2, b_idx2);
				faces.resize(4);
				faces[0].set(s.a, s.b, s.c, true);
				faces[1].set(s.b, s.c, s.d, false);
				faces[2].set(s.a, s.b, s.d, false);
				faces[3].set(s.a, s.c, s.d, false);
			}
			break;
		case 3://triangle contains origin
			{
				Vector3f n=(s.b.p-s.a.p).cross(s.c.p-s.a.p);
				Vector3f point=support(A, n, a_idx)-support(B, -n, b_idx);
				Vector3f point2=support(A, -n, a_idx2)-support(B, n, b_idx2);
				if(a_idx==s.a.a_idx&&b_idx==s.a.b_idx||a_idx==s.b.a_idx&&b_idx==s.b.b_idx||a_idx==s.c.a_idx&&b_idx==s.c.b_idx)
				{
				//	float u, v, w;
					depth=point.magnitude();
					normal=n;
					float mag=normal.mag_sq();
					if(mag)
					{
						normal*=1/sqrt(mag);
						get_uvw_face(s.a.p, s.b.p, s.c.p, u, v, w);
						ca=u*A.w[s.a.a_idx]+v*A.w[s.b.a_idx]+w*A.w[s.c.a_idx];
						cb=u*B.w[s.a.b_idx]+v*B.w[s.b.b_idx]+w*B.w[s.c.b_idx];
					}
					else//abc are collinear
					{
						normal=-point;
						normal*=1/depth;
						//normal=s.b.p-s.a.p;
						//normal*=1/normal.magnitude();
						//make_orthogonal(normal);
						get_uv_line(s.a.p, s.b.p, u, v);
						ca=u*A.w[s.a.a_idx]+v*A.w[s.b.a_idx];
						cb=u*B.w[s.a.b_idx]+v*B.w[s.b.b_idx];
					}
					return true;
				}
				if(a_idx2==s.a.a_idx&&b_idx2==s.a.b_idx||a_idx2==s.b.a_idx&&b_idx2==s.b.b_idx||a_idx2==s.c.a_idx&&b_idx2==s.c.b_idx)
				{
				//	float u, v, w;
					depth=point2.magnitude();
					normal=n;
					float mag=normal.mag_sq();
					if(mag)
					{
						normal*=1/sqrt(mag);
						get_uvw_face(s.a.p, s.b.p, s.c.p, u, v, w);
						ca=u*A.w[s.a.a_idx]+v*A.w[s.b.a_idx]+w*A.w[s.c.a_idx];
						cb=u*B.w[s.a.b_idx]+v*B.w[s.b.b_idx]+w*B.w[s.c.b_idx];
					}
					else//abc are collinear
					{
						normal=-point2;
						normal*=1/depth;
						//normal=s.b.p-s.a.p;
						//normal*=1/normal.magnitude();
						//make_orthogonal(normal);
						get_uv_line(s.a.p, s.b.p, u, v);
						ca=u*A.w[s.a.a_idx]+v*A.w[s.b.a_idx];
						cb=u*B.w[s.a.b_idx]+v*B.w[s.b.b_idx];
					}
					return true;
				}
				s.add(point, a_idx, b_idx);
				PolytopeVertex p5(point2, a_idx2, b_idx2);
				faces.resize(6);
				faces[0].set(s.a, s.b, s.d, false);
				faces[1].set(s.a, s.c, s.d, false);
				faces[2].set(s.b, s.c, s.d, false);
				faces[3].set(s.a, s.b, p5, false);
				faces[4].set(s.a, s.c, p5, false);
				faces[5].set(s.b, s.c, p5, false);
			}
			break;
		case 4:
			faces.resize(4);
			faces[0].set(s.a, s.b, s.c, false);
			faces[1].set(s.b, s.c, s.d, false);
			faces[2].set(s.a, s.b, s.d, false);
			faces[3].set(s.a, s.c, s.d, false);
			break;
		}
		//switch(s.count)//expand the simplex to tetrahedron
		//{
		//case 1://point
		//	{
		//		Vector3f point=support(A, s.a.p, a_idx)-support(B, -s.a.p, b_idx);
		//		s.add(point, a_idx, b_idx);
		//	}//continued in line case
		//case 2://line segment
		//	{
		//		Vector3f edge=s.b.p-s.a.p, n=edge.triple_product(s.b.p, edge);
		//		Vector3f point=support(A, n, a_idx)-support(B, -n, b_idx);
		//		if(a_idx==s.a.a_idx&&b_idx==s.a.b_idx||a_idx==s.b.a_idx&&b_idx==s.b.b_idx)//farthest point is already in simplex
		//		{
		//			depth=point.magnitude();
		//			//normal=s.b.p-s.a.p;
		//			//normal*=1/normal.magnitude();
		//			//make_orthogonal(normal);
		//			normal=-point;
		//			normal*=1/depth;

		//			get_uv_line(s.a.p, s.b.p, u, v);
		//			ca=u*A.w[s.a.a_idx]+v*A.w[s.b.a_idx];
		//			cb=u*B.w[s.a.b_idx]+v*B.w[s.b.b_idx];
		//			return true;
		//		}
		//		s.add(point, a_idx, b_idx);
		//	}//continued in triangle case
		//case 3://triangle
		//	{
		//		Vector3f n=(s.b.p-s.a.p).cross(s.c.p-s.a.p);
		//		Vector3f point=support(A, n, a_idx)-support(B, -n, b_idx);
		//		if(a_idx==s.a.a_idx&&b_idx==s.a.b_idx||a_idx==s.b.a_idx&&b_idx==s.b.b_idx||a_idx==s.c.a_idx&&b_idx==s.c.b_idx)
		//		{
		//		//	float u, v, w;
		//			depth=point.magnitude();
		//			normal=(s.b.p-s.a.p).cross(s.c.p-s.a.p);
		//			float mag=normal.mag_sq();
		//			if(mag)
		//			{
		//				normal*=1/sqrt(mag);
		//				get_uvw_face(s.a.p, s.b.p, s.c.p, u, v, w);
		//				ca=u*A.w[s.a.a_idx]+v*A.w[s.b.a_idx]+w*A.w[s.c.a_idx];
		//				cb=u*B.w[s.a.b_idx]+v*B.w[s.b.b_idx]+w*B.w[s.c.b_idx];
		//			}
		//			else//abc are collinear
		//			{
		//				normal=-point;
		//				normal*=1/depth;
		//				//normal=s.b.p-s.a.p;
		//				//normal*=1/normal.magnitude();
		//				//make_orthogonal(normal);
		//				get_uv_line(s.a.p, s.b.p, u, v);
		//				ca=u*A.w[s.a.a_idx]+v*A.w[s.b.a_idx];
		//				cb=u*B.w[s.a.b_idx]+v*B.w[s.b.b_idx];
		//			}
		//			return true;
		//		}
		//		s.add(point, a_idx, b_idx);
		//	}
		//	break;
		//}

		//std::vector<Face> faces(4);
		//faces[0].set(s.a, s.b, s.c, s.count<=3);
		//faces[1].set(s.b, s.c, s.d, false);
		//faces[2].set(s.a, s.b, s.d, false);
		//faces[3].set(s.a, s.c, s.d, false);
		int boundary_face=0;
		unsigned final_face=0;
		float final_distance;
		if(A.idx==g_a_idx&&B.idx==g_b_idx)
			g_faces.push_back(faces), ++n_EPA_iterations;//

		for(int iteration=0;iteration<64;++iteration)
		{
			//for(unsigned k=0;k<faces.size();++k)
			//	if(faces[k].n.isnan_or_inf())
			//		int LOL_1=0;
			final_face=0, final_distance=faces[0].origin_distance;
			for(int k=1, kEnd=faces.size();k<kEnd;++k)//find the closest face to the origin
				if(final_distance>faces[k].origin_distance)
					final_distance=faces[k].origin_distance, final_face=k;
			auto &f=faces[final_face];

			Vector3f point=support(A, f.n, a_idx)-support(B, -f.n, b_idx);//get farthest point in that direction
			//Vector3f point2;
			//int a_idx2=-1, b_idx2=-1;
			//bool contains_origin=f.contains_origin;
			//if(contains_origin)//get farthest point in both directions
			//	point2==support(A, -f.n, a_idx2)-support(B, f.n, b_idx2);
			
			//for(int k=0, kEnd=faces.size();k<kEnd;++k)//check if new point already here
			//{
			//	auto &f2=faces[k];
			//	if(f2.v1.a_idx==a_idx&&f2.v1.b_idx==b_idx||f2.v2.a_idx==a_idx&&f2.v2.b_idx==b_idx||f2.v3.a_idx==a_idx&&f2.v3.b_idx==b_idx)
			//		int LOL_1=0;
			//}
			//if(f.v1.a_idx==a_idx&&f.v1.b_idx==b_idx||f.v2.a_idx==a_idx&&f.v2.b_idx==b_idx||f.v3.a_idx==a_idx&&f.v3.b_idx==b_idx)//if it is part of this face then terminate
			//	break;
			bool boundary=false;
			for(int k=0, kEnd=faces.size();k<kEnd;++k)//check if new point already here
			{
				auto &f2=faces[k];
				if(f2.v1.a_idx==a_idx&&f2.v1.b_idx==b_idx||f2.v2.a_idx==a_idx&&f2.v2.b_idx==b_idx||f2.v3.a_idx==a_idx&&f2.v3.b_idx==b_idx)
				{
				//	boundary=true;
					boundary_face=final_face, boundary=true;
				//	boundary_face=k, boundary=true;
					break;
				}
			}
			if(boundary)
				break;
			float D=point.dot(f.n);
			if(D-f.origin_distance<0.00001)
			{
				boundary_face=final_face;
				break;
			}
			//	break;

			//insert the new point
			typedef std::pair<PolytopeVertex, PolytopeVertex> PolytopeEdge;
			std::vector<PolytopeEdge> loose_edges;
			for(unsigned k=0;k<faces.size();)//remove faces not facing the new point
			{
				auto &f=faces[k];
				if(f.n.dot(point-f.v1.p)>0)//triangle k faces the new point, remove it
				{
					for(int ke=0;ke<3;++ke)//add this face's edges to the loose edge list
					{
						auto &v1=(&f.v1)[ke], &v2=(&f.v1)[(ke+1)%3];
						bool new_edge=true;
						for(unsigned ke2=0;ke2<loose_edges.size();++ke2)//check if the edge is already there 
						{
							auto &v3=loose_edges[ke2].first;
							auto &v4=loose_edges[ke2].second;
							if(v1.a_idx==v3.a_idx&&v1.b_idx==v3.b_idx&&v2.a_idx==v4.a_idx&&v2.b_idx==v4.b_idx
								||v1.a_idx==v4.a_idx&&v1.b_idx==v4.b_idx&&v2.a_idx==v3.a_idx&&v2.b_idx==v3.b_idx)
							{
								loose_edges.erase(loose_edges.begin()+ke2);
								new_edge=false;
								break;
							}
						}
						if(new_edge)
							loose_edges.push_back(PolytopeEdge(v1, v2));
					}
					faces.erase(faces.begin()+k), final_face-=final_face>k;
				}
				else
					++k;
			}
			PolytopeVertex v(point, a_idx, b_idx);
			int n_faces=faces.size(), n_edges=loose_edges.size();
			faces.resize(n_faces+n_edges);
			for(int k=0;k<n_edges;++k)//reconstruct polytope
				faces[n_faces+k].set(v, loose_edges[k].first, loose_edges[k].second, false);
			//faces.resize(faces.size()+3);
			//auto &f2=faces[final_face];
			//int s=faces.size()-3;
			//PolytopeVertex v(point, a_idx, b_idx);
			//faces[s].set(v, f2.v2, f2.v3);
			//faces[s+1].set(f2.v1, v, f2.v3);
			//faces[s+2].set(f2.v1, f2.v2, v);
			//faces.erase(faces.begin()+final_face);//make sure polytope is convex
			if(A.idx==g_a_idx&&B.idx==g_b_idx)
				g_faces.push_back(faces), ++n_EPA_iterations;//
		}

	/*	std::vector<PolytopeVertex> pt(4);//polytope
		pt[0]=s.a, pt[1]=s.b, pt[2]=s.c, pt[3]=s.d;
		int index;
		for(;;)
		{
			const float TOLERANCE=1e-10f;
			float min_distance=infinity;
			Vector3f n2;

			for(int k=0, kEnd=pt.size();k<kEnd;++k)//get nearest edge
			{
				int k2=(k+1)%kEnd;
				Vector3f &v1=pt[k].p, &v2=pt[k2].p;
				Vector3f edge=v2-v1;//for each edge in polytope
				Vector3f n=v1*(edge.dot(edge))-edge*(edge.dot(v1));//(AxB)xC = B(C.dot(A))?A(C.dot(B))	triple product to get vector from edge to the origin
				n/=n.magnitude();
				float dist=n.dot(v1);
				if(min_distance>dist)
				{
					min_distance=dist;
					index=k, n2=n;
				}
			}
			Vector3f sup=support(A, n2, a_idx)-support(B, -n2, b_idx);
			float d=sup.dot(n2);
			if(d-min_distance<=TOLERANCE)
			{
				normal=n2;
				depth=min_distance;
				break;
			}
			else//insert sup
			{
				bool duplicate=false;
				for(int k=0, kEnd=pt.size();k<kEnd;++k)
					if(pt[k].p==sup)
					{
						duplicate=true;
						break;
					}
				if(duplicate)//
				{
					normal=n2;
					depth=min_distance;
					not_sure=true;
					break;
				}
				pt.insert(pt.begin()+(index+1)%pt.size(), PolytopeVertex(sup, a_idx, b_idx));
			}
		}//*/
		//find contact point
		//if(frame_number==13)
		//	int LOL_1=0;
	//	g_faces=faces;//
		auto &f=faces[boundary_face];
		//float mag_n=f.n.magnitude();
		//depth=f.origin_distance/mag_n;
		//normal=f.n/mag_n;
		if(f.contains_origin)
		{
			depth=0;
			if(s0_count==1)//point
			{
				normal.set(0, 0, 1);
				ca=A.w[s.a.a_idx];
				cb=B.w[s.a.b_idx];
			}
			else//line
			{
				normal=s.b.p-s.a.p;
				normal*=1/normal.magnitude();
				make_orthogonal(normal);
				get_uv_line(s.a.p, s.b.p, u, v);
				ca=u*A.w[s.a.a_idx]+v*A.w[s.b.a_idx];
				cb=u*B.w[s.a.b_idx]+v*B.w[s.b.b_idx];
			}
		}
		else
		{
			depth=f.origin_distance;
		//	normal=(f.v2.p-f.v1.p).cross(f.v3.p-f.v1.p);
		//	normal*=1/normal.magnitude();
			normal=f.n;
		//	float u, v, w;
			get_uvw_face(f.v1.p, f.v2.p, f.v3.p, u, v, w);
			vector_to_origin=u*f.v1.p+v*f.v2.p+w*f.v3.p;
		//	m->draw_line(Vector3f(0, 0, 0), u*f.v1.p+v*f.v2.p+w*f.v3.p, 0x0000FF);//
			ca=u*A.w[f.v1.a_idx]+v*A.w[f.v2.a_idx]+w*A.w[f.v3.a_idx];
			cb=u*B.w[f.v1.b_idx]+v*B.w[f.v2.b_idx]+w*B.w[f.v3.b_idx];
		}
	/*	int index_v2=(index+1)%pt.size();
		PolytopeVertex &v1=pt[index], &v2=pt[index_v2];
		int ka1=v1.a_idx, ka2=v2.a_idx,
			kb1=v1.b_idx, kb2=v2.b_idx;
		float u, v;
		Vector3f n, p;
		get_uv_line(v1.p, v2.p, u, v, n);
		if(u<=0)
			ca=A.w[v2.a_idx], cb=B.w[v2.b_idx];//, p=v2.p;
		else if(v<=0)
			ca=A.w[v1.a_idx], cb=B.w[v1.b_idx];//, p=v1.p;
		else
			ca=u*A.w[v1.a_idx]+v*A.w[v2.a_idx], cb=u*B.w[v1.b_idx]+v*B.w[v2.b_idx];//, p=u*v1.p+v*v2.p;//*/
		return true;
	}
	return false;
}

void			find_least_penetration(PhysObject const &A, PhysObject const &B, double &max_dist, int &Aface_idx, int &Bv_idx, Vector3f &normal)
{
	for(unsigned kf=0;kf<A.tr.size();++kf)
	{
		auto &tr=A.tr[kf];
		Vector3f const &v1=A.w[tr.a], &v2=A.w[tr.b], &v3=A.w[tr.c];
		Vector3f n=-(v2-v1).cross(v3-v1);//points into A
		n/=n.magnitude();
		int kv;
		Vector3f vb=support(B, n, kv);
		double dist=(vb-v1).dot(n);
		if(max_dist<dist)
			max_dist=dist, Aface_idx=kf, Bv_idx=kv, normal=-n;//points out of A
	}
}
//bool			is_minkowski_face(Vector3f const &a, Vector3f const &b, Vector3f const &c, Vector3f const &d, Vector3f const &b_x_a, Vector3f const d_x_c)//normals of adjacend faces of the edges
////bool			is_minkowski_face(Vector3f const &a, Vector3f const &b, Vector3f const &c, Vector3f const &d)
//{
////	Vector3f b_x_a=b.cross(a), d_x_c=d.cross(c);//test is arcs ab and cd intersect on the unit sphere
//	float cba=c.dot(b_x_a),
//		dba=d.dot(b_x_a),
//		adc=a.dot(d_x_c),
//		bdc=b.dot(d_x_c);
//	return cba*dba<0&&adc*bdc<0&&cba*bdc<0;
//}
bool			build_minkowski_face(PhysObject const &A, int ea, PhysObject const &B, int eb)
{
	auto &edge_a=A.e[ea], &edge_b=B.e[eb];
//	auto &a_tr1=A.tr[edge_a.tr1], &a_tr2=A.tr[edge_a.tr2], &b_tr1=B.w[edge_b.tr1], &b_tr2=B.w[edge_b.tr2];
	Vector3f a=A.get_normal(edge_a.tr1),
		b=A.get_normal(edge_a.tr2),
		c=-B.get_normal(edge_b.tr1),
		d=-B.get_normal(edge_b.tr2);

	//is_minkowski_face - test is arcs ab and cd intersect on the unit sphere
	Vector3f b_x_a=A.w[edge_a.v2]-A.w[edge_a.v1], d_x_c=B.w[edge_b.v2]-B.w[edge_b.v1];
	//Vector3f b_x_a=b.cross(a), d_x_c=d.cross(c);
	float cba=c.dot(b_x_a),
		dba=d.dot(b_x_a),
		adc=a.dot(d_x_c),
		bdc=b.dot(d_x_c);
	return cba*dba<0&&adc*bdc<0&&cba*bdc<0;
}
float			SAT_distance(PhysObject const &A, int ea, PhysObject const &B, int eb)
{
	auto &edge_a=A.e[ea], &edge_b=A.e[eb];
	const Vector3f &av1=A.w[edge_a.v1], &bv1=B.w[edge_b.v1];
	Vector3f adir=A.w[edge_a.v2]-av1, bdir=B.w[edge_b.v2]-bv1;
	if(adir.cross(bdir).mag_sq()==0)//skip parallel edges
		return -infinity;
	Vector3f normal=adir.cross(bdir);
	normal*=1/normal.magnitude();
	if(normal.dot(av1-A.p)>0)//assume normal points A->B
		normal=-normal;
	return normal.dot(bv1-av1);
}
int				manifold_support(std::vector<Vector3f> const &out, Vector3f const &direction)
{
	int index=0;
	double max_dot=out[0].dot(direction);
	for(int k=1, kEnd=out.size();k<kEnd;++k)
	{
		double dot=out[k].dot(direction);
		if(max_dot<dot)
			max_dot=dot, index=k;
	}
	return index;
}
//struct			Edge
//{
//	Point p0, p1;
//	Edge():p0(0, 0), p1(0, 0){}
//	Edge(Point const &p0, Point const &p1):p0(p0), p1(p1){}
//	void set(Point const &p0, Point const &p1){this->p0=p0, this->p1=p1;}
//};
//bool			intersection_SAT(PhysPolygon const &A, PhysPolygon const &B, double &depth, Point &normal, Point *ca, int &count_a, Point *cb, int &count_b)
bool			intersection_SAT(PhysObject const &A, PhysObject const &B, double &depth, Vector3f &normal, Vector3f *cp, int &count)
{
	double Amax_dist=-infinity;
	int Aface_idx, Bv_idx;
	Vector3f an;
	find_least_penetration(A, B, Amax_dist, Aface_idx, Bv_idx, an);
	if(Amax_dist>0)
		return false;
	double Bmax_dist=-infinity;
	int Bface_idx, Av_idx;
	Vector3f bn;
	find_least_penetration(B, A, Bmax_dist, Bface_idx, Av_idx, bn);
	if(Bmax_dist>0)
		return false;

	int edge_a=-1, edge_b=-1;
	float edge_dist=-infinity;
	for(int ka=0, kaEnd=A.e.size();ka<kaEnd;++ka)//query edge direction
	{
		auto &ea=A.e[ka];
		for(int kb=0, kbEnd=B.e.size();kb<kbEnd;++kb)
		{
			auto &eb=B.e[kb];
			if(build_minkowski_face(A, ka, B, kb))
			{
				float separation=SAT_distance(A, ka, B, kb);
				if(edge_dist<separation)
					edge_dist=separation, edge_a=ka, edge_b=kb;
			}
		}
	}

//	count=0;


/*	bool swapped=Amax_dist<Bmax_dist;//A's side: largest negative penetration
	PhysPolygon const *A2, *B2;
	if(swapped)
	{
		A2=&B, B2=&A;
		Point const &b0=B.w[Bface_idx], &b1=B.w[(Bface_idx+1)%B.w.size()];
		depth=-Bmax_dist, normal.set(b0.y-b1.y, b1.x-b0.x);//pointing out of A
	//	draw_line(b0, b0+normal);
	}
	else
	{
		A2=&A, B2=&B;
		Point const &a0=A.w[Aface_idx], &a1=A.w[(Aface_idx+1)%A.w.size()];
		depth=-Amax_dist, normal.set(a1.y-a0.y, a0.x-a1.x);//pointing out of A
	//	draw_line(a0, a0+normal);
	}
	normal/=normal.magnitude();
	//Edge reference(A2->w[Aface_idx], A2->w[(Aface_idx+1)%A2->w.size()]);
	Edge incident (B2->w[Bface_idx], B2->w[(Bface_idx+1)%B2->w.size()]);
	//Point const &a0=A2->w[Aface_idx], &a1=A2->w[(Aface_idx+1)%A2->w.size()];
	//Point const &b0=B2->w[Bface_idx], &b1=B2->w[(Bface_idx+1)%B2->w.size()];
	//reference.set(a0, a1);
	//incident.set(b0, b1);
	bool incident_exists=true;
	Edge in(incident);
//	Edge in(incident), out;
	if(frame_number==14&&A.idx==1&&B.idx==5)
		int LOL_1=0;
	for(int ka=0, a_size=A2->w.size();ka<a_size;++ka)
	{
		Edge a(A2->w[ka], A2->w[(ka+1)%a_size]);
		double d_b0=distance_to_line_times_seg_length(a.p1, a.p0, in.p0),
			d_b1=distance_to_line_times_seg_length(a.p1, a.p0, in.p1);
		if(d_b0>0)
		{
			if(d_b1>0)//both inside, save b1
				;
			else//b0 inside, b1 outside, save i
			{
				Point i;
				intersection_lines_nonparallel(a.p0, a.p1, in.p0, in.p1, i);
				in.p1=i;
			}
		}
		else if(d_b1>0)//b0 outside, b1 inside, store i & b1
		{
			Point i;
			intersection_lines_nonparallel(a.p0, a.p1, in.p0, in.p1, i);
			in.p1=i;
		}
		else//both outside, none saved
		{
			incident_exists=false;
			break;
		}
		//if(d_b0>0)
		//{
		//	if(d_b1>0)//both inside, save b1
		//	{
		//		out.p0=in.p0;
		//		out.p1=in.p1;
		//	}
		//	else//b0 inside, b1 outside, save i
		//	{
		//		Point i;
		//		intersection_lines_nonparallel(a0, a1, b0, b1, i);
		//		out.p0=in.p0;
		//		out.p1=i;
		//	}
		//}
		//else if(d_b1>0)//b0 outside, b1 inside, store i & b1
		//{
		//	Point i;
		//	intersection_lines_nonparallel(a0, a1, b0, b1, i);
		//	out.p0=i;
		//	out.p1=in.p1;
		//}
		////else both outside, none saved
		//in=out;
	}
	if(incident_exists)
	{
		count=2, cp[0]=in.p0, cp[1]=in.p1;
		draw_line(cp[0]-10, cp[0]+10), draw_line(cp[0].x-10, cp[0].y+10, cp[0].x+10, cp[0].y-10);//
		draw_line(cp[1]-10, cp[1]+10), draw_line(cp[1].x-10, cp[1].y+10, cp[1].x+10, cp[1].y-10);//
	}
	else
	{
		count=0;
		GUIPrint(ghMemDC, A.p.x, A.p.y, "no incident");
		GUIPrint(ghMemDC, B.p.x, B.p.y, "no incident");
		draw_line(incident.p0+10, incident.p1+10);
	}//*/

/*	bool As_side=Amax_dist>Bmax_dist;//largest negative penetration
	int reference, incident;
	if(As_side)
	{
		reference=Aface_idx, incident=Bface_idx;
		Point const &a0=A.w[Aface_idx], &a1=A.w[(Aface_idx+1)%A.w.size()];
		depth=-Amax_dist, normal.set(a1.y-a0.y, a0.x-a1.x);//pointing out of A
	//	draw_line(a0, a0+normal);
	}
	else
	{
		reference=Bface_idx, incident=Aface_idx;
		Point const &b0=B.w[Bface_idx], &b1=B.w[(Bface_idx+1)%B.w.size()];
		depth=-Bmax_dist, normal.set(b0.y-b1.y, b1.x-b0.x);//pointing out of A
	//	draw_line(b0, b0+normal);
	}
	normal/=normal.magnitude();

	std::vector<Point> in=B.w, out;//crop of B inside A - Sutherland?Hodgman algorithm
	for(int ka=0, a_size=A.w.size();ka<a_size;++ka)
	{
		Point const &a0=A.w[ka], &a1=A.w[(ka+1)%a_size];
		for(unsigned kb=0, b_size=in.size();kb<b_size;++kb)
		{
			Point const &b0=in[kb], &b1=in[(kb+1)%b_size];
			double d_b0=distance_to_line_times_seg_length(a1, a0, b0),
				d_b1=distance_to_line_times_seg_length(a1, a0, b1);
			if(d_b0>0)
			{
				if(d_b1>0)//both inside, save b1
					out.push_back(b1);
				else//b0 inside, b1 outside, save i
				{
					Point i;
					intersection_lines_nonparallel(a0, a1, b0, b1, i);
					out.push_back(i);
				}
			}
			else if(d_b1>0)//b0 outside, b1 inside, store i & b1
			{
				Point i;
				intersection_lines_nonparallel(a0, a1, b0, b1, i);
				out.push_back(i);
				out.push_back(b1);
			}
			//else both outside, none saved
		}
		in=std::move(out);
		//HPEN hPen=CreatePen(PS_DOT, 1, 0xFF0000);
		//hPen=(HPEN)SelectObject(ghMemDC, hPen);
		//for(int ki=0, i_size=in.size();ki<i_size;++ki)
		//{
		//	Point p0=2*(in[ki]-B.p)+w/2-ka*100,
		//		p1=2*(in[(ki+1)%i_size]-B.p)+w/2-ka*100;
		//	draw_line(p0, p1);
		//	GUIPrint(ghMemDC, p0.x, p0.y, ki);
		//}
		//hPen=(HPEN)SelectObject(ghMemDC, hPen);
		//DeleteObject(hPen);
	}
	//GUIPrint(ghMemDC, B.idx*75., h/2+100.+A.idx*18, "%d vs %d", A.idx, B.idx);
	//HPEN hPen=CreatePen(PS_DOT, 1, 0xFF0000);
	//hPen=(HPEN)SelectObject(ghMemDC, hPen);
	//for(int ki=0, i_size=in.size();ki<i_size;++ki)
	//{
	//	Point p0=4*(in[ki]-B.p), p1=4*(in[(ki+1)%i_size]-B.p);
	//	draw_line(p0.x+B.idx*75, p0.y+h/2, p1.x+B.idx*75, p1.y+h/2);
	//	GUIPrint(ghMemDC, p0.x+B.idx*75, p0.y+h/2, ki);
	////	Point p0=in[ki]+2, p1=in[(ki+1)%i_size]+2;
	//	//draw_line(p0, p1);
	//	//GUIPrint(ghMemDC, p0.x, p0.y, ki);
	//}
	//hPen=(HPEN)SelectObject(ghMemDC, hPen);
	//DeleteObject(hPen);

	if(in.size()==1)
		cp[0]=in[0], count=1;
	else if(in.size()>1)
	{
		Point const *a0, *a1;
		if(As_side)
			a0=&A.w[Aface_idx], a1=&A.w[(Aface_idx+1)%A.w.size()];
		else
			a0=&B.w[Bface_idx], a1=&B.w[(Bface_idx+1)%B.w.size()];
		count=2;
		cp[0]=in[manifold_support(in, *a0-*a1)];
		cp[1]=in[manifold_support(in, *a1-*a0)];
		//draw_line(cp[0]-10, cp[0]+10), draw_line(cp[0].x-10, cp[0].y+10, cp[0].x+10, cp[0].y-10);
		//draw_line(cp[1]-10, cp[1]+10), draw_line(cp[1].x-10, cp[1].y+10, cp[1].x+10, cp[1].y-10);
	}
	else count=0;//*/

	//if(info)
	//{
	//	HPEN hPen=CreatePen(PS_DOT, 1, 0xFF0000);
	//	hPen=(HPEN)SelectObject(ghMemDC, hPen);
	//	GUIPrint(ghMemDC, A.p.x, A.p.y, "intersect");
	//	GUIPrint(ghMemDC, B.p.x, B.p.y, "intersect");
	//	for(int kp=0;kp<count;++kp)
	//		draw_line(A.p, cp[kp]);
	//	//draw_line(B.w[Bv_idx], B.w[Bv_idx]+an);
	//	//draw_line(A.w[Av_idx], A.w[Av_idx]+bn);
	//	hPen=(HPEN)SelectObject(ghMemDC, hPen);
	//	DeleteObject(hPen);
	//	draw_line(cp[0], cp[0]+normal*100);
	//}
	return true;
}

int				find_contact(PhysObject const &A, PhysObject const &B)
{
	for(int kc=0, kcEnd=contacts.size();kc<kcEnd;++kc)//find contact with these objects
	{
		auto &c=contacts[kc];
		if(c.a_idx==A.idx&&c.b_idx==B.idx||c.a_idx==B.idx&&c.b_idx==A.idx)
			return kc;
	}
	return -1;
}
void			insert_contact_point(PhysObject const &A, PhysObject const &B, Vector3f const &n, Vector3f const &t1, Vector3f const &t2, Vector3f const &ra, Vector3f const &rb, float depth, int kc, bool &new_contact)
{
	bool persistent=false;
	int kp=0;
	if(warm_starting&&kc!=-1)
	{
		auto &c=contacts[kc];
		for(int kpEnd=c.p.size();kp<kpEnd;++kp)//search for close previous contact
		{
			auto &p=c.p[kp];
			Vector3f dra=ra-p.ra, drb=rb-p.rb;
			if(abs(dra.x)<proximity_limit&&abs(dra.y)<proximity_limit&&abs(drb.x)<proximity_limit&&abs(drb.y)<proximity_limit)
			{
				persistent=true;
				break;
			}
		}
	}
	if(persistent)
	{
		auto &c=contacts[kc];
	//	c.persistent=true;//
		c.p[kp].old=false;
		c.p[kp].ra=ra;
		c.p[kp].rb=rb;
		c.p[kp].depth=depth;
	}
	else if(kc!=-1)
		contacts[kc].add_cp(ra, rb, depth);
	else if(new_contact)
	{
		contacts.push_back(ContactInfo(A.idx, B.idx, depth, ra, rb, n, t1, t2));
		new_contact=false;
	}
	else
		contacts.rbegin()->add_cp(ra, rb, depth);
}
int				find_contact(PhysObject const &A, Vector3f const &n)
{
	for(int kc=0, kcEnd=contacts.size();kc<kcEnd;++kc)//find A's contact with boundary
	{
		auto &c=contacts[kc];
		if(c.a_idx==A.idx&&c.b_idx==-2&&c.n==n)
			return kc;
	}
	return -1;
}
void			insert_contact_point(PhysObject const &A, Vector3f const &n, Vector3f const &t1, Vector3f const &t2, Vector3f const &ra, Vector3f const &rb, float depth, int kc, bool &new_contact)
{
	bool persistent=false;
	int kp=0;
	if(warm_starting&&kc!=-1)
	{
		auto &c=contacts[kc];
		for(int kpEnd=c.p.size();kp<kpEnd;++kp)
		{
			auto &p=c.p[kp];
			Vector3f dra=ra-p.ra,//local
				drb=rb-p.rb;//world, point may have moved on the surface
			if(abs(dra.x)<proximity_limit&&abs(dra.y)<proximity_limit&&abs(drb.dot(n))<proximity_limit)
			{
				persistent=true;
				break;
			}
		}
	}
	if(persistent)
	{
		auto &c=contacts[kc];
	//	c.persistent=true;//
		c.p[kp].old=false;
		c.p[kp].ra=ra;
		c.p[kp].rb=rb;
		c.p[kp].depth=depth;
	}
	else if(kc!=-1)
		contacts[kc].add_cp(ra, rb, depth);
	else if(new_contact)
	{
		contacts.push_back(ContactInfo(A.idx, -2, depth, ra, rb, n, t1, t2));//normal pointing out of A
		new_contact=false;
	}
	else
		contacts.rbegin()->add_cp(ra, rb, depth);
}
struct			ImpulseInfo
{
	Vector3f pa, pb,//position
		i,//impulse
		r;//rotational impulse
	bool boundary;
//	ImpulseInfo(){}
//	ImpulseInfo(Vector3f const &pa, Vector3f const &i, Vector3f const &r):boundary(true), pa(pa), i(i), r(r){}
//	ImpulseInfo(Vector3f const &pa, Vector3f const &pb, Vector3f const &i, Vector3f const &r):boundary(false), pa(pa), pb(pb), i(i), r(r){}
};
//std::vector<std::vector<ImpulseInfo>> impulse_info;//
std::vector<ImpulseInfo> impulse_info;//
std::vector<ContactInfo> contacts0;//
void			render()
{
	if(timer){		// if(kb['Q'])						for(int k=0;k<npols;++k){auto &p=env[k];if(p.tx_idx>=0)p.tx1.x+=rand()%3-1, p.tx1.y+=rand()%3-1, p.tx2.x+=rand()%3-1, p.tx2.y+=rand()%3-1, p.tx3.x+=rand()%3-1, p.tx3.y+=rand()%3-1, transform(p);}
					 if(kb[VK_SHIFT]){	 if(kb['W'])	cam.x+=10*dcam*cax*cay,	cam.y+=10*dcam*sax*cay,	cam.z+=10*dcam*say;
										 if(kb['A'])	cam.x-=10*dcam*sax,		cam.y+=10*dcam*cax;
										 if(kb['S'])	cam.x-=10*dcam*cax*cay,	cam.y-=10*dcam*sax*cay,	cam.z-=10*dcam*say;
										 if(kb['D'])	cam.x+=10*dcam*sax,		cam.y-=10*dcam*cax;
									//	 if(kb['T'])	cam.x-=10.*dcam*cax*say,	cam.y-=10.*dcam*sax*say,	cam.z+=10.*dcam*cay;
									//	 if(kb['G'])	cam.x+=10.*dcam*cax*say,	cam.y+=10.*dcam*sax*say,	cam.z-=10.*dcam*cay;}
										 if(kb['T'])	cam.z+=10*dcam;
										 if(kb['G'])	cam.z-=10*dcam;}
				else{					 if(kb['W'])	cam.x+=dcam*cax*cay,		cam.y+=dcam*sax*cay,		cam.z+=dcam*say;
										 if(kb['A'])	cam.x-=dcam*sax,			cam.y+=dcam*cax;
										 if(kb['S'])	cam.x-=dcam*cax*cay,		cam.y-=dcam*sax*cay,		cam.z-=dcam*say;
										 if(kb['D'])	cam.x+=dcam*sax,			cam.y-=dcam*cax;
									//	 if(kb['T'])	cam.x-=dcam*cax*say,		cam.y-=dcam*sax*say,		cam.z+=dcam*cay;
									//	 if(kb['G'])	cam.x+=dcam*cax*say,		cam.y+=dcam*sax*say,		cam.z-=dcam*cay;}
										 if(kb['T'])	cam.z+=dcam;
										 if(kb['G'])	cam.z-=dcam;}
					 if(kb[VK_UP]){						ay+=da_tfov*da, cay=cos(ay), say=sin(ay);	for(;ay<0;)ay+=_2pi;	for(;ay>_2pi;)ay-=_2pi;}
					 if(kb[VK_DOWN]){					ay-=da_tfov*da, cay=cos(ay), say=sin(ay);	for(;ay<0;)ay+=_2pi;	for(;ay>_2pi;)ay-=_2pi;}
					 if(kb[VK_LEFT]){					ax+=da_tfov*da, cax=cos(ax), sax=sin(ax);	for(;ax<0;)ax+=_2pi;	for(;ax>_2pi;)ax-=_2pi;}
					 if(kb[VK_RIGHT]){					ax-=da_tfov*da, cax=cos(ax), sax=sin(ax);	for(;ax<0;)ax+=_2pi;	for(;ax>_2pi;)ax-=_2pi;}
					 if(kb[VK_ADD])						tanfov/=DTANFOV, fov/=DFOV, m->changeFov(), da_tfov=tanfov>1?1:tanfov;
					 if(kb[VK_SUBTRACT])				tanfov*=DTANFOV, fov*=DFOV, m->changeFov(), da_tfov=tanfov>1?1:tanfov;
					 if(kb[VK_RBUTTON]){				int boom=10;
														for(int k=0;k<npols;++k){	auto &p=env[k];
																					p.p1.x+=(rand()%boom)-boom/2, p.p1.y+=(rand()%boom)-boom/2, p.p1.z+=(rand()%boom)-boom/2;
																					p.p2.x+=(rand()%boom)-boom/2, p.p2.y+=(rand()%boom)-boom/2, p.p2.z+=(rand()%boom)-boom/2;
																					p.p3.x+=(rand()%boom)-boom/2, p.p3.y+=(rand()%boom)-boom/2, p.p3.z+=(rand()%boom)-boom/2, transform(p);}}}
	m->clear_screen();

	int gcount=10;
	float gstep=100, gz=0.1f;
	//float gstep=abs(cam.z)+10;
	float gsize=gcount*gstep, px=gstep*floor(cam.x/gstep), py=gstep*floor(cam.y/gstep);//draw floor grid
	for(int kx=-gcount;kx<=gcount;++kx)
	{
		float xg=px+gstep*kx;
		m->draw_line(Vector3f(xg, py-gsize, gz), Vector3f(xg, py+gsize, gz), 0);
	}
	for(int ky=-gcount;ky<=gcount;++ky)
	{
		float yg=py+gstep*ky;
		m->draw_line(Vector3f(px-gsize, yg, gz), Vector3f(px+gsize, yg, gz), 0);
	}
	if(tick)
	{
		for(int k=0, kEnd=objects.size();k<kEnd;++k)//simulate
		{
			auto &A=objects[k];
			A.iterate(timescale);
		}
		//if(frame_number==276)
		//if(frame_number==37)
		//if(frame_number==32)
		if(frame_number==574)
			int LOL_1=0;
		for(int ka=0, kaEnd=objects.size()-1;ka<kaEnd;++ka)//collision detection
		{
			auto &A=objects[ka];
			for(int kb=ka+1, kbEnd=objects.size();kb<kbEnd;++kb)
			{
				auto &B=objects[kb];
				float r=A.r+B.r;
				if(abs(A.p.x-B.p.x)<=r&&abs(A.p.y-B.p.y)<=r&&abs(A.p.z-B.p.z)<=r)
				{
					Vector3f normal, ca, cb;
					float depth;
					if(intersection_EPA(A, B, normal, depth, ca, cb))
				//	Simplex s;
				//	if(intersection_polyhedra_GJK(A, B, A.p-B.p, s))
					{
						g_ca=ca, g_cb=cb;//
						if(info)
						{
							m->print(ca, "%f", depth);
							m->print(cb, "%f", depth);
							m->draw_line(A.p, ca, 0);
							m->draw_line(B.p, cb, 0);
						//	draw_projected_yz(A, 0), draw_projected_yz(B, 0);
						//	draw_line_3projections(A.p, ca, 0), draw_line_3projections(B.p, cb, 0);
						//	draw_line_projected_yz(A.p, ca, 0), draw_line_projected_yz(B.p, cb, 0);
							//m->print(A.p, "intersect");
							//m->print(B.p, "intersect");
						}
						Vector3f ra=A.world_to_local(ca), rb=B.world_to_local(cb);
						Vector3f t1(normal.y, -normal.x, 0), t2=normal.cross(t1);
						int kc=find_contact(A, B);
						bool new_contact=true;
						insert_contact_point(A, B, normal, t1, t2, ra, rb, depth, kc, new_contact);
					//	contacts.push_back(ContactInfo(ka, kb, depth, ra, rb, normal, t1, t2));
						if(A.idx==g_a_idx&&B.idx==g_b_idx)
							g_c=kc==-1?contacts.size()-1:kc;
					}
				}
			}
		}
		for(int k=0, kEnd=objects.size();k<kEnd;++k)//detect contact with boundaries
		{
			auto &A=objects[k];
			//m->print(A.p, "%d", A.idx);//
			if(A.p.z-A.r<zground)//ground
			{
				Vector3f n_ground(0, 0, -1), t1(1, 0, 0), t2(0, 1, 0);
				int kc=find_contact(A, n_ground);
				bool new_contact=kc==-1;
				for(int kv=0, kvEnd=A.w.size();kv<kvEnd;++kv)
				{
					if(A.w[kv].z<zground)
					{
						float depth=zground-A.w[kv].z;
						Vector3f ra=A.world_to_local(A.w[kv]), rb(A.w[kv].x, A.w[kv].y, zground);
						insert_contact_point(A, n_ground, t1, t2, ra, rb, depth, kc, new_contact);
					}
				}
			}
			if(A.p.x-A.r<xleft)//x=0 wall (yz plane)
			{
				Vector3f n(-1, 0, 0), t1(0, 1, 0), t2(0, 0, 1);
				int kc=find_contact(A, n);
				bool new_contact=kc==-1;
				for(int kv=0, kvEnd=A.w.size();kv<kvEnd;++kv)
				{
					if(A.w[kv].x<xleft)
					{
						float depth=xleft-A.w[kv].x;
						Vector3f ra=A.world_to_local(A.w[kv]), rb(xleft, A.w[kv].y, A.w[kv].z);
						insert_contact_point(A, n, t1, t2, ra, rb, depth, kc, new_contact);
					}
				}
			}
			if(A.p.y-A.r<yfront)//y=0 wall (xz plane)
			{
				Vector3f n(0, -1, 0), t1(1, 0, 0), t2(0, 0, 1);
				int kc=find_contact(A, n);
				bool new_contact=kc==-1;
				for(int kv=0, kvEnd=A.w.size();kv<kvEnd;++kv)
				{
					if(A.w[kv].y<yfront)
					{
						float depth=yfront-A.w[kv].y;
						Vector3f ra=A.world_to_local(A.w[kv]), rb(A.w[kv].x, yfront, A.w[kv].z);
						insert_contact_point(A, n, t1, t2, ra, rb, depth, kc, new_contact);
					}
				}
			}
		}
		if(frame_number==33)
			int LOL_1=0;
		for(unsigned kc=0;kc<contacts.size();)//remove old contacts
		{
			auto &c=contacts[kc];
			for(unsigned kp=0;kp<c.p.size();)
			{
				auto &p=c.p[kp];
				if(c.b_idx==-2)//with boundary
				{
					if(p.old)
						c.p.erase(c.p.begin()+kp);
					else
						++kp;
				}
				else if(c.b_idx>=0)
				{
					auto &A=objects[c.a_idx], &B=objects[c.b_idx];
					Vector3f pa=A.local_to_world(p.ra), pb=B.local_to_world(p.rb);
					float depth=(pa-pb).dot(c.n);
					if(depth<=0||!inclusion_point_convex_polyhedron(A, pb)&&!inclusion_point_convex_polyhedron(B, pa))//check if pa is inside B and pb is inside A
						c.p.erase(c.p.begin()+kp);
					else
						++kp;
				}
			}
			if(c.p.size())
				++kc;
			else
				contacts.erase(contacts.begin()+kc);
		}
		for(int k=0, kEnd=objects.size();k<kEnd;++k)//apply forces
		{
			auto &A=objects[k];
			A.v.z-=gravity*timescale;
		}

		contacts0=contacts;//
		for(int kc=0, kcEnd=contacts.size();kc<kcEnd;++kc)//resolve contacts: pre step
		{
			auto &c=contacts[kc];
			for(int kp=0, kpEnd=c.p.size();kp<kpEnd;++kp)
			{
				auto &p=c.p[kp];
				p.old=true;
				if(c.b_idx==-2)//boundary
				{
					auto &A=objects[c.a_idx];
					Vector3f &n=c.n, &t1=c.t1, &t2=c.t2;
					Vector3f ra=A.local_to_world_rotate_only(p.ra);
					float ra_n=ra.dot(n), ra_t1=ra.dot(t1), ra_t2=ra.dot(t2), ra2=ra.dot(ra);
					
					p.mass_n =1/(A.inv_mass+(A.inv_I*ra.cross(n )).cross(ra).dot(n ));
					p.mass_t1=1/(A.inv_mass+(A.inv_I*ra.cross(t1)).cross(ra).dot(t1));
					p.mass_t2=1/(A.inv_mass+(A.inv_I*ra.cross(t2)).cross(ra).dot(t2));

					//Vector3f ra_x_n=ra.cross(n), ra_x_t1=ra.cross(t1), ra_x_t2=ra.cross(t2);//v4
					//p.mass_n =1/(A.inv_mass+(A.inv_I*ra_x_n ).dot(ra_x_n ));
					//p.mass_t1=1/(A.inv_mass+(A.inv_I*ra_x_t1).dot(ra_x_t1));
					//p.mass_t2=1/(A.inv_mass+(A.inv_I*ra_x_t2).dot(ra_x_t2));

					//p.mass_n =1/(A.inv_mass+(ra2-ra_n *ra_n )*A.inv_I.determinant_diagonal());//v5
					//p.mass_t1=1/(A.inv_mass+(ra2-ra_t1*ra_t1)*A.inv_I.determinant_diagonal());//
					//p.mass_t2=1/(A.inv_mass+(ra2-ra_t2*ra_t2)*A.inv_I.determinant_diagonal());//

					p.bias=-beta/timescale*clamp_negative(-p.depth+penetration_slop);
					//m->print(0, h/2+18*point, "%.2f, %.2f, %.2f, %.2f", p.mass_n, p.mass_t1, p.mass_t2, p.bias);//
					if(accumulated_impulses)
					{
						Vector3f P=p.Pn*n+p.Pt1*t1+p.Pt2*t2;
						A.v-=A.inv_mass*P,		A.vr-=A.inv_I*ra.cross(P);
					}
				}
				else if(c.b_idx>=0)
				{
					auto &A=objects[c.a_idx], &B=objects[c.b_idx];
					Vector3f &n=c.n, &t1=c.t1, &t2=c.t2;
					Vector3f ra=A.local_to_world_rotate_only(p.ra), rb=B.local_to_world_rotate_only(p.rb);
					float ra_n=ra.dot(n), ra_t1=ra.dot(t1), ra_t2=ra.dot(t2), rb_n=rb.dot(n), rb_t1=rb.dot(t1), rb_t2=rb.dot(t2),
						ra2=ra.dot(ra), rb2=rb.dot(rb);
					float ma_mb=A.inv_mass+B.inv_mass;
					
					p.mass_n =1/(ma_mb+((A.inv_I*ra.cross(n )).cross(ra)+(B.inv_I*rb.cross(n )).cross(rb)).dot(n ));//https://pybullet.org/Bullet/phpBB3/viewtopic.php?t=8886
					p.mass_t1=1/(ma_mb+((A.inv_I*ra.cross(t1)).cross(ra)+(B.inv_I*rb.cross(t1)).cross(rb)).dot(t1));
					p.mass_t2=1/(ma_mb+((A.inv_I*ra.cross(t2)).cross(ra)+(B.inv_I*rb.cross(t2)).cross(rb)).dot(t2));
					//Vector3f
					//	ra_x_n=ra.cross(n), ra_x_t1=ra.cross(t1), ra_x_t2=ra.cross(t2),
					//	rb_x_n=rb.cross(n), rb_x_t1=rb.cross(t1), rb_x_t2=rb.cross(t2);
					//p.mass_n =1/(ma_mb+(A.inv_I*ra_x_n ).dot(ra_x_n )+(B.inv_I*rb_x_n ).dot(rb_x_n ));//v4
					//p.mass_t1=1/(ma_mb+(A.inv_I*ra_x_t1).dot(ra_x_t1)+(B.inv_I*rb_x_t1).dot(rb_x_t1));
					//p.mass_t1=1/(ma_mb+(A.inv_I*ra_x_t2).dot(ra_x_t2)+(B.inv_I*rb_x_t2).dot(rb_x_t2));

					//p.mass_n =1/(ma_mb+(ra2-ra_n *ra_n )*A.inv_I.determinant_diagonal()+(rb2-rb_n *rb_n )*B.inv_I.determinant_diagonal());//v5
					//p.mass_t1=1/(ma_mb+(ra2-ra_t1*ra_t1)*A.inv_I.determinant_diagonal()+(rb2-rb_t1*rb_t1)*B.inv_I.determinant_diagonal());
					//p.mass_t2=1/(ma_mb+(ra2-ra_t2*ra_t2)*A.inv_I.determinant_diagonal()+(rb2-rb_t2*rb_t2)*B.inv_I.determinant_diagonal());

					p.bias=-beta/timescale*clamp_negative(-p.depth+penetration_slop);
					//m->print(0, h/2+18*point, "%.2f, %.2f, %.2f, %.2f", p.mass_n, p.mass_t1, p.mass_t2, p.bias);//
					if(accumulated_impulses)
					{
						Vector3f P=p.Pn*n+p.Pt1*t1+p.Pt2*t2;
						A.v-=A.inv_mass*P;		A.vr-=A.inv_I*ra.cross(P);
						A.v+=B.inv_mass*P;		B.vr+=B.inv_I*rb.cross(P);
					}
				}
			}
		}
		int ncp=0;
		impulse_info.resize(0);
		for(unsigned kc=0;kc<contacts.size();++kc)
			ncp+=contacts[kc].p.size();
		impulse_info.resize(ncp);//call constructor for all
		//if(frame_number==31)
		if(frame_number==575)//0 vs 6
			int LOL_1=0;
		for(int i=0;i<n_iterations;++i)//apply impulses
		{
			int cp=0;
			for(int kc=0, kcEnd=contacts.size();kc<kcEnd;++kc)
			{
				auto &c=contacts[kc];
				for(int kp=0, kpEnd=c.p.size();kp<kpEnd;++kp, ++cp)
				{
					auto &p=c.p[kp];
					if(c.b_idx==-2)//with boundary
					{
						auto &A=objects[c.a_idx];
						Vector3f &n=c.n, &t1=c.t1, &t2=c.t2;
						Vector3f ra=A.local_to_world_rotate_only(p.ra);
						Vector3f dv=-A.v-A.vr.cross(ra);
						float vn=dv.dot(n);
						float dPn=p.mass_n*(-vn+p.bias);
						if(accumulated_impulses)
						{
							float Pn0=p.Pn;
							p.Pn=clamp_positive(Pn0+dPn);
							dPn=p.Pn-Pn0;
						}
						else
							dPn=clamp_positive(dPn);
						Vector3f Pn=dPn*n;
						A.v-=Pn*A.inv_mass;			A.vr-=A.inv_I*ra.cross(Pn);
								
						float f_limit=sqrt(boundary_friction*A.friction)*dPn;

						dv=-A.v-A.vr.cross(ra);
						float vt1=dv.dot(t1);
						float dPt1=p.mass_t1*-vt1;
						if(accumulated_impulses)
						{
							float oldTangentImpulse=p.Pt1;
							p.Pt1=clamp(oldTangentImpulse+dPt1, -f_limit, f_limit);
							dPt1=p.Pt1-oldTangentImpulse;
						}
						else
							dPt1=clamp(dPt1, -f_limit, f_limit);
						Vector3f Pt1=dPt1*t1;
						A.v-=Pt1*A.inv_mass;		A.vr-=A.inv_I*ra.cross(Pt1);

						dv=-A.v-A.vr.cross(ra);
						float vt2=dv.dot(t2);
						float dPt2=p.mass_t2*-vt2;
						if(accumulated_impulses)
						{
							float oldTangentImpulse=p.Pt2;
							p.Pt2=clamp(oldTangentImpulse+dPt2, -f_limit, f_limit);
							dPt2=p.Pt2-oldTangentImpulse;
						}
						else
							dPt2=clamp(dPt2, -f_limit, f_limit);
						Vector3f Pt2=dPt2*t2;
						A.v-=Pt2*A.inv_mass;		A.vr-=A.inv_I*ra.cross(Pt2);
					//	A.v-=(Pt1+Pt2)*A.inv_mass,	A.vr-=A.inv_I*ra.cross(Pt1+Pt2);

						auto &ii=impulse_info[cp];
						ii.boundary=true, ii.pa=A.local_to_world(p.ra), ii.i+=Pn+Pt1+Pt2;
						if(info)
						{
							Vector3f cp=A.local_to_world(p.ra);
							m->draw_line(cp, cp-(Pn+Pt1+Pt2), 0);//
						}
					}
					else if(c.b_idx>=0)//with object
					{
						auto &A=objects[c.a_idx], &B=objects[c.b_idx];
						Vector3f &n=c.n, &t1=c.t1, &t2=c.t2;
						Vector3f ra=A.local_to_world_rotate_only(p.ra), rb=B.local_to_world_rotate_only(p.rb);
						Vector3f dv=B.v+B.vr.cross(rb)-A.v-A.vr.cross(ra);
						float vn=dv.dot(n);
						float dPn=p.mass_n*(-vn+p.bias);
						if(accumulated_impulses)
						{
							float Pn0=p.Pn;
							p.Pn=clamp_positive(Pn0+dPn);
							dPn=p.Pn-Pn0;
						}
						else
							dPn=clamp_positive(dPn);
						Vector3f Pn=dPn*n;
						A.v-=Pn*A.inv_mass;			A.vr-=A.inv_I*ra.cross(Pn);
						B.v+=Pn*B.inv_mass;			B.vr+=B.inv_I*rb.cross(Pn);
						
						float f_limit=sqrt(A.friction*B.friction)*dPn;

						dv=B.v+B.vr.cross(rb)-A.v-A.vr.cross(ra);
						float vt1=dv.dot(t1);
						float dPt1=-vt1*p.mass_t1;
						if(accumulated_impulses)
						{
							float oldTangentImpulse=p.Pt1;
							p.Pt1=clamp(oldTangentImpulse+dPt1, -f_limit, f_limit);
							dPt1=p.Pt1-oldTangentImpulse;
						}
						else
							dPt1=clamp(dPt1, -f_limit, f_limit);
						Vector3f Pt1=dPt1*t1;
						A.v-=Pt1*A.inv_mass;		A.vr-=A.inv_I*ra.cross(Pt1);
						B.v+=Pt1*B.inv_mass;		B.vr+=B.inv_I*rb.cross(Pt1);

						dv=B.v+B.vr.cross(rb)-A.v-A.vr.cross(ra);
						float vt2=dv.dot(t2);
						float dPt2=-vt2*p.mass_t2;
						if(accumulated_impulses)
						{
							float oldTangentImpulse=p.Pt2;
							p.Pt2=clamp(oldTangentImpulse+dPt2, -f_limit, f_limit);
							dPt2=p.Pt2-oldTangentImpulse;
						}
						else
							dPt2=clamp(dPt2, -f_limit, f_limit);
						Vector3f Pt2=dPt2*t2;
						A.v-=Pt2*A.inv_mass;		A.vr-=A.inv_I*ra.cross(Pt2);
						B.v+=Pt2*B.inv_mass;		B.vr+=B.inv_I*rb.cross(Pt2);

						auto &ii=impulse_info[cp];
						ii.boundary=false, ii.pa=A.local_to_world(p.ra), ii.pb=B.local_to_world(p.rb), ii.i+=Pn+Pt1+Pt2;
						if(info)
						{
							Vector3f cpa=A.local_to_world(p.ra), cpb=A.local_to_world(p.rb);//
							m->draw_line(cpa, cpa-(Pn+Pt1+Pt2), 0), m->draw_line(cpb, cpb+(Pn+Pt1+Pt2), 0);//
						}
					}
				}
			}
		}
#if 0
		for(int kc=0, kcEnd=contacts.size();kc<kcEnd;++kc)//resolve contacts
		{
			auto &c=contacts[kc];
			if(c.b_idx==-2)//boundary
			{
				auto &A=objects[c.a_idx];
				float C_F=sqrt(A.friction*boundary_friction);
				for(int kp=0, kpEnd=c.p.size();kp<kpEnd;++kp)
				{
					auto &p=c.p[kp];
					p.old=true;
					Vector3f ra=A.local_to_world_rotate_only(p.ra);
					Vector3f raxn=ra.cross(c.n), raxt1=ra.cross(c.t1), raxt2=ra.cross(c.t2);
					Vector3f invM_JnT [2]={-c.n *A.inv_mass, A.I_diagonal?A.inv_I.multiply_diagonal(-raxn ):A.inv_I*(-raxn )};
					Vector3f invM_Jt1T[2]={-c.t1*A.inv_mass, A.I_diagonal?A.inv_I.multiply_diagonal(-raxt1):A.inv_I*(-raxt1)};
					Vector3f invM_Jt2T[2]={-c.t2*A.inv_mass, A.I_diagonal?A.inv_I.multiply_diagonal(-raxt2):A.inv_I*(-raxt2)};
					float eff_mass_n=-c.n.dot(invM_JnT[0])-raxn.dot(invM_JnT[1]);
					float eff_mass_t1=-c.t1.dot(invM_Jt1T[0])-raxt1.dot(invM_Jt1T[1]);
					float eff_mass_t2=-c.t2.dot(invM_Jt2T[0])-raxt2.dot(invM_Jt2T[1]);

					float normalImpulseSum=0, t1ImpulseSum=0, t2ImpulseSum=0;

					float Jn_V1=-c.n.dot(A.v)-raxn.dot(A.vr);
					float Jt1_V1=-c.t1.dot(A.v)-raxt1.dot(A.vr);
					float Jt2_V1=-c.t2.dot(A.v)-raxt2.dot(A.vr);
					float bn=-beta/timescale*p.depth+C_R*(-A.v-A.vr.cross(ra)).dot(c.n);
					float bt1=-beta/timescale*p.depth;
					float bt2=-beta/timescale*p.depth;
					float lambda_n=-(Jn_V1+bn)/eff_mass_n;
					float lambda_t1=-(Jt1_V1+bt1)/eff_mass_t1;
					float lambda_t2=-(Jt2_V1+bt2)/eff_mass_t2;

					float Pn=normalImpulseSum;//clamp normal impulse
					normalImpulseSum=clamp_positive(normalImpulseSum+lambda_n);
					lambda_n=normalImpulseSum-Pn;
					
					float Pt1=t1ImpulseSum;//clamp t1 impulse
					float t_limit=C_F*lambda_n;
					t1ImpulseSum=clamp(t1ImpulseSum+lambda_t1, -t_limit, t_limit);
					lambda_t1=t1ImpulseSum-Pt1;

					float Pt2=t2ImpulseSum;//clamp t2 impulse
					t2ImpulseSum=clamp(t2ImpulseSum+lambda_t2, -t_limit, t_limit);
					lambda_t2=t2ImpulseSum-Pt2;

					A.v +=invM_JnT[0]*lambda_n+invM_Jt1T[0]*lambda_t1+invM_Jt2T[0]*lambda_t2;
					A.vr+=invM_JnT[1]*lambda_n+invM_Jt1T[1]*lambda_t1+invM_Jt2T[1]*lambda_t2;
					impulse_info[A.idx].push_back(ImpulseInfo(A.local_to_world(p.ra), invM_JnT[0]*lambda_n+invM_Jt1T[0]*lambda_t1+invM_Jt2T[0]*lambda_t2,
						invM_JnT[1]*lambda_n+invM_Jt1T[1]*lambda_t1+invM_Jt2T[1]*lambda_t2));
				}
			}
			else if(c.b_idx>=0)
			{
				auto &A=objects[c.a_idx], &B=objects[c.b_idx];
				float C_F=sqrt(A.friction*B.friction);
				for(int kp=0, kpEnd=c.p.size();kp<kpEnd;++kp)
				{
					auto &p=c.p[kp];
					p.old=true;
					Vector3f ra=A.local_to_world_rotate_only(p.ra), rb=B.local_to_world_rotate_only(p.rb);
					Vector3f raxn=ra.cross(c.n), raxt1=ra.cross(c.t1), raxt2=ra.cross(c.t2);
					Vector3f rbxn=rb.cross(c.n), rbxt1=rb.cross(c.t1), rbxt2=rb.cross(c.t2);
					Vector3f invM_JnT [4]={-c.n *A.inv_mass, A.I_diagonal?A.inv_I.multiply_diagonal(-raxn ):A.inv_I*(-raxn ), c.n *B.inv_mass, B.I_diagonal?B.inv_I.multiply_diagonal(rbxn ):B.inv_I*rbxn };
					Vector3f invM_Jt1T[4]={-c.t1*A.inv_mass, A.I_diagonal?A.inv_I.multiply_diagonal(-raxt1):A.inv_I*(-raxt1), c.t1*B.inv_mass, B.I_diagonal?B.inv_I.multiply_diagonal(rbxt1):B.inv_I*rbxt1};
					Vector3f invM_Jt2T[4]={-c.t2*A.inv_mass, A.I_diagonal?A.inv_I.multiply_diagonal(-raxt2):A.inv_I*(-raxt2), c.t2*B.inv_mass, B.I_diagonal?B.inv_I.multiply_diagonal(rbxt2):B.inv_I*rbxt2};
					float eff_mass_n =-c.n .dot(invM_JnT [0])-raxn .dot(invM_JnT [1])+c.n .dot(invM_JnT [2])+rbxn .dot(invM_JnT [3]);
					float eff_mass_t1=-c.t1.dot(invM_Jt1T[0])-raxt1.dot(invM_Jt1T[1])+c.t1.dot(invM_Jt1T[2])+rbxt1.dot(invM_Jt1T[3]);
					float eff_mass_t2=-c.t2.dot(invM_Jt2T[0])-raxt2.dot(invM_Jt2T[1])+c.t2.dot(invM_Jt2T[2])+rbxt2.dot(invM_Jt2T[3]);
					
					float normalImpulseSum=0, t1ImpulseSum=0, t2ImpulseSum=0;

					float Jn_V1 =-c.n .dot(A.v)-raxn .dot(A.vr)+c.n .dot(B.v)+rbxn .dot(B.vr);
					float Jt1_V1=-c.t1.dot(A.v)-raxt1.dot(A.vr)+c.t1.dot(B.v)+rbxt1.dot(B.vr);
					float Jt2_V1=-c.t2.dot(A.v)-raxt2.dot(A.vr)+c.t2.dot(B.v)+rbxt2.dot(B.vr);
					float bn=-beta/timescale*p.depth+C_R*(-A.v-A.vr.cross(ra)+B.v+B.vr.cross(rb)).dot(c.n);
					float bt1=-beta/timescale*p.depth;
					float bt2=-beta/timescale*p.depth;
					float lambda_n=-(Jn_V1+bn)/eff_mass_n;
					float lambda_t1=-(Jt1_V1+bt1)/eff_mass_t1;
					float lambda_t2=-(Jt2_V1+bt2)/eff_mass_t2;

					float Pn=normalImpulseSum;//clamp normal impulse
					normalImpulseSum=clamp_positive(normalImpulseSum+lambda_n);
					lambda_n=normalImpulseSum-Pn;
					
					float Pt1=t1ImpulseSum;//clamp t1 impulse
					float t_limit=C_F*lambda_n;
					t1ImpulseSum=clamp(t1ImpulseSum+lambda_t1, -t_limit, t_limit);
					lambda_t1=t1ImpulseSum-Pt1;

					float Pt2=t2ImpulseSum;//clamp t2 impulse
					t2ImpulseSum=clamp(t2ImpulseSum+lambda_t2, -t_limit, t_limit);
					lambda_t2=t2ImpulseSum-Pt2;

					A.v +=invM_JnT[0]*lambda_n+invM_Jt1T[0]*lambda_t1+invM_Jt2T[0]*lambda_t2;
					A.vr+=invM_JnT[1]*lambda_n+invM_Jt1T[1]*lambda_t1+invM_Jt2T[1]*lambda_t2;
					B.v +=invM_JnT[2]*lambda_n+invM_Jt1T[2]*lambda_t1+invM_Jt2T[2]*lambda_t2;
					B.vr+=invM_JnT[3]*lambda_n+invM_Jt1T[3]*lambda_t1+invM_Jt2T[3]*lambda_t2;
					impulse_info[A.idx].push_back(ImpulseInfo(A.local_to_world(p.ra), invM_JnT[0]*lambda_n+invM_Jt1T[0]*lambda_t1+invM_Jt2T[0]*lambda_t2,
						invM_JnT[1]*lambda_n+invM_Jt1T[1]*lambda_t1+invM_Jt2T[1]*lambda_t2));//
					impulse_info[B.idx].push_back(ImpulseInfo(B.local_to_world(p.rb), invM_JnT[2]*lambda_n+invM_Jt1T[2]*lambda_t1+invM_Jt2T[2]*lambda_t2,
						invM_JnT[3]*lambda_n+invM_Jt1T[3]*lambda_t1+invM_Jt2T[3]*lambda_t2));//
				}
			}
		}
#endif
		if(!warm_starting)
			contacts.resize(0);
		++frame_number;
	}
	float total_energy=0;
	if(info)
	{
		for(int ko=0, koEnd=objects.size();ko<koEnd;++ko)//calculate total energy
		{
			auto &A=objects[ko];
			float kinetic_energy=0.5f*(A.mass*A.v.mag_sq() + A.vr.dot(A.I*A.vr)), potential_energy=A.mass*gravity*(A.p.z-zground);
			total_energy+=kinetic_energy+potential_energy;
		}
		for(unsigned ki=0;ki<impulse_info.size();++ki)//draw impulses
		{
			auto &ii=impulse_info[ki];
			m->draw_line(ii.pa, ii.pa-ii.i, 0);
			if(!ii.boundary)
				m->draw_line(ii.pb, ii.pb+ii.i, 0);
		}
		if(g_c!=-1&&g_c<(int)contacts0.size())//draw normals
		{
			auto &c=contacts0[g_c];
			if(c.a_idx==g_a_idx&&c.b_idx==g_b_idx)
			{
				auto &A=objects[c.a_idx], &B=objects[c.b_idx];
				for(int kp=0, kpEnd=c.p.size();kp<kpEnd;++kp)
				{
					auto &p=c.p[kp];
					Vector3f cpa=A.local_to_world(p.ra);
					m->draw_line(cpa, cpa+c.n*10, 0);
				}
			}
		}
		for(int kc=0, kcEnd=contacts0.size();kc<kcEnd;++kc)//draw contacts
		{
			auto &c=contacts0[kc];
			for(int kp=0, kpEnd=c.p.size();kp<kpEnd;++kp)
			{
				auto &p=c.p[kp];
				if(c.b_idx==-2)//with boundary
				{
					auto &A=objects[c.a_idx];
					Vector3f ra=A.local_to_world(p.ra);
					m->draw_line(A.p, ra, 0);
					m->draw_line(ra, ra+p.Pn+p.Pt1+p.Pt2, 0);
					if(p.old)
						m->print(ra, "old");
				}
				else if(c.b_idx>=0)//with object
				{
					auto &A=objects[c.a_idx], &B=objects[c.b_idx];
					Vector3f ra=A.local_to_world(p.ra), rb=B.local_to_world(p.rb);
					m->draw_line(A.p, ra, 0), m->draw_line(B.p, rb, 0);
					m->draw_line(ra, ra+p.Pn+p.Pt1+p.Pt2, 0), m->draw_line(rb, rb-(p.Pn+p.Pt1+p.Pt2), 0);
					if(p.old)
						m->print(ra, "old"), m->print(rb, "old");
				}
			}
		}
	}
	Vector3f zero(0, 0, 0);
//	m->draw_line(0, 0, (float)w, (float)h, 0);
	m->draw_line(zero, Vector3f(10000, 0, 0), 0);//draw axes
	m->draw_line(zero, Vector3f(0, 10000, 0), 0);
	m->draw_line(zero, Vector3f(0, 0, 10000), 0);

//	m->draw_ground();
	m->render();

	if(info)
	{
		//{
		//	draw_line_3projections(zero, Vector3f(100, 0, 0), 0);//draw projections
		//	draw_line_3projections(zero, Vector3f(0, 100, 0), 0);
		//	draw_line_3projections(zero, Vector3f(0, 0, 100), 0);
		//	//draw_line_projected_yz(zero, Vector3f(0, 1000, 0), 0);
		//	//draw_line_projected_yz(zero, Vector3f(0, 0, 1000), 0);
		//	for(int k=0, kEnd=objects.size();k<kEnd;++k)//
		//		draw_3projections(objects[k], 0);
		//	//	draw_projected_yz(objects[k], 0);//
		//}
		for(int k=0, kEnd=objects.size();k<kEnd;++k)
			m->print(objects[k].p, "%d", objects[k].idx);
		if(objects.size()>=2)
		{
			draw_minkowski_difference(objects[g_a_idx], objects[g_b_idx], 0);
			m->draw_line(Vector3f(0, 0, 0), vector_to_origin, 0x0000FF);
			if(EPA_iteration>=0&&EPA_iteration<g_faces.size())
			{
				for(int k=0, kEnd=g_faces[EPA_iteration].size();k<kEnd;++k)
				{
					auto &f=g_faces[EPA_iteration][k];
					char intensity=255*k/kEnd;
					char color[4]={255-intensity, intensity, 0, 0};
					m->draw_triangle(f.v1.p, f.v2.p, f.v3.p, *(int*)color, true);
					//m->draw_triangle(f.v1.p, f.v2.p, f.v3.p, rand()<<15|rand(), true);
					//m->draw_line(f.v1.p, f.v2.p, 0x00FF00);
					//m->draw_line(f.v2.p, f.v3.p, 0x00FF00);
					//m->draw_line(f.v3.p, f.v1.p, 0x00FF00);
				}
			}
			//Vector3f zero(0, 0, 0);
			//if(g_ca!=zero&&g_cb!=zero)
			//{
			//	m->draw_line(objects[0].p, g_ca, 0x00FF00);
			//	m->draw_line(objects[1].p, g_cb, 0x0000FF);//rgb
			//}
			m->print(0, 18, "%d of %d GJK iterations, %d EPA iterations", EPA_iteration, n_GJK_iterations, n_EPA_iterations-n_GJK_iterations);
		}
	}
	m->print(0, 18*2, "4,5: new scenario {unpaused,paused}");
	m->print(0, 18*3, "R,Q: reset {unpaused,paused}");
	m->print(0, 18*4, "E/ENTER: next frame");
	m->print(0, 18*5, "P,F: toggle pause");
	m->print(0, 18*6, "Y: toggle info");
	m->print(0, 18*7, "timescale: %f", timescale);
	m->print(0, 18*8, "message: 0x%08x", g_message);
	if(info)
		m->print(0, 18*9, "Total energy: %f", total_energy);
	//if(impulse_info.size()==objects.size())
	//{
	//	for(unsigned k=0;k<objects.size();++k)
	//	{
	//		auto &io=impulse_info[k];
	//		for(unsigned ki=0;ki<io.size();++ki)
	//		{
	//			m->draw_line(io[ki].p, io[ki].p+100*io[ki].i, 0xFF0000);
	//			m->draw_line(io[ki].p, io[ki].p+100*io[ki].r, 0x00FF00);
	//		}
	//	}
	//}
	m->print(w*2/5, 0, "frame #%08d, scenario #%03d", frame_number-tick, scenario_number);
//	GUIPrint(ghDC, w*2/5, 0, "frame #%08d, scenario #%03d", frame_number-1, scenario_number);
	tick=false;
	m->show();
}