#include <windows.h>
#include <gl\GL.h>
#include <gl\GLU.h>
#include <gl\glut.h>
#include <iostream>

#include "semiLag.hpp"
#include "semiLag.cpp"
#include "vector.hpp"
#include "vector.cpp"
//#define particle_Grid_Len 10
//
//const double h = 1.0 / particle_Grid_Len;
//
///* 粒子结构 */
//
//struct Particle
//{
//	unsigned int  r, g, b;      /* 粒子的颜色 */
//	double vx, vy, vz;          /* 粒子的当前速度 */
//	double duration;			/* 粒子持续时间 */
//	double size;                /* 粒子尺寸 */
//	double temperature;         /* 粒子温度 */
//	double density;             /* 粒子密度 */
//	double force[3];
//};
//
///* 粒子网格类 */
//class Grid
//{
//public:
//	Particle***   particle;               /* 粒子指针 */
//	//int         particle_Cnt;         /* 粒子数目 */
//	double      temp_avg;            /* 温度条件 */
//
//	/* 构造函数 */
//	Grid() {
//		particle = NULL;
//		//particle_Cnt = 0;
//		temp_avg = 0;
//	}
//	/* 析构函数 */
//	~Grid() {
//		for (int i = 0; i < particle_Grid_Len; i++)
//			for (int j = 0; j < particle_Grid_Len; j++) {
//				delete[] particle[i][j];
//				particle[i][j] = NULL;
//			}
//
//		for (int i = 0; i < particle_Grid_Len; i++) {
//			delete[] particle[i];
//			particle[i] = NULL;
//		}
//		delete[] particle;
//		particle = NULL;
//	}
//
//	/* 创建粒子数组 */
//	void CreateGrid(long num) {
//		if (particle)
//			delete[] particle;
//
//		particle = (Particle***)new Particle**[particle_Grid_Len];
//
//		for (int i = 0; i < particle_Grid_Len; i++) {
//			particle[i] = (Particle**)new Particle*[particle_Grid_Len];
//		}
//
//		for (int i = 0; i < particle_Grid_Len; i++)
//			for (int j = 0; j < particle_Grid_Len; j++) {
//				particle[i][j] = new Particle[particle_Grid_Len];
//			}
//
//
//	}
//
//	/* 设置和获取颜色属性 */
//	void SetAllColor(GLint r, GLint g, GLint b) {
//		for (int i = 0; i < particle_Grid_Len; i++)
//			for (int j = 0; j < particle_Grid_Len; j++)
//				for (int k = 0; k < particle_Grid_Len; k++) {
//					particle[i][j][k].r = r;
//					particle[i][j][k].g = g;
//					particle[i][j][k].b = b;
//				}
//	}
//
//	bool SetColor(GLint i, GLint j, GLint k, GLint r, GLint g, GLint b) {
//		if (i >= 0 && i < particle_Grid_Len
//			&& j >= 0 && j < particle_Grid_Len
//			&& k >= 0 && k < particle_Grid_Len) {
//			particle[i][j][k].r = r;
//			particle[i][j][k].g = g;
//			particle[i][j][k].b = b;
//			return true;
//		}
//		return false;
//	}
//
//	bool GetColor(GLint i, GLint j, GLint k, GLint &r, GLint &g, GLint &b) {
//		if (i >= 0 && i < particle_Grid_Len
//			&& j >= 0 && j < particle_Grid_Len
//			&& k >= 0 && k < particle_Grid_Len) {
//			r = particle[i][j][k].r;
//			g = particle[i][j][k].g;
//			b = particle[i][j][k].b;
//			return true;
//		}
//		return false;
//	}
//
//	/* 设置和获取速度属性 */
//	void SetAllVelocity(GLdouble vx, GLdouble vy, GLdouble vz) {
//		for (int i = 0; i < particle_Grid_Len; i++)
//			for (int j = 0; j < particle_Grid_Len; j++)
//				for (int k = 0; k < particle_Grid_Len; k++) {
//					particle[i][j][k].vx = vx;
//					particle[i][j][k].vy = vy;
//					particle[i][j][k].vz = vz;
//				}
//	}
//
//	bool SetVelocity(GLint i, GLint j, GLint k, GLdouble vx, GLdouble vy, GLdouble vz) {
//		if (i >= 0 && i < particle_Grid_Len
//			&& j >= 0 && j < particle_Grid_Len
//			&& k >= 0 && k < particle_Grid_Len) {
//			particle[i][j][k].vx = vx;
//			particle[i][j][k].vy = vy;
//			particle[i][j][k].vz = vz;
//			return true;
//		}
//		return false;
//	}
//
//	bool GetVelocity(GLint i, GLint j, GLint k, GLdouble &vx, GLdouble &vy, GLdouble &vz) {
//		if (i >= 0 && i < particle_Grid_Len
//			&& j >= 0 && j < particle_Grid_Len
//			&& k >= 0 && k < particle_Grid_Len) {
//			vx = particle[i][j][k].vx;
//			vy = particle[i][j][k].vy;
//			vz = particle[i][j][k].vz;
//			return true;
//		}
//		return false;
//	}
//
//	/* 设置和获取持续时间 */
//	void SetAllDuration(GLdouble duration) {
//		for (int i = 0; i < particle_Grid_Len; i++)
//			for (int j = 0; j < particle_Grid_Len; j++)
//				for (int k = 0; k < particle_Grid_Len; k++) {
//					particle[i][j][k].duration = duration;
//				}
//	}
//
//	bool SetDuration(GLint i, GLint j, GLint k, GLdouble duration) {
//		if (i >= 0 && i < particle_Grid_Len
//			&& j >= 0 && j < particle_Grid_Len
//			&& k >= 0 && k < particle_Grid_Len) {
//			particle[i][j][k].duration = duration;
//			return true;
//		}
//		return false;
//	}
//
//	bool GetDuration(GLint i, GLint j, GLint k, GLdouble &duration) {
//		if (i >= 0 && i < particle_Grid_Len
//			&& j >= 0 && j < particle_Grid_Len
//			&& k >= 0 && k < particle_Grid_Len) {
//			duration = particle[i][j][k].duration;
//			return true;
//		}
//		return false;
//	}
//
//	/* 设置和获取尺寸属性 */
//	void SetAllSize(GLdouble size) {
//		for (int i = 0; i < particle_Grid_Len; i++)
//			for (int j = 0; j < particle_Grid_Len; j++)
//				for (int k = 0; k < particle_Grid_Len; k++) {
//					particle[i][j][k].size = size;
//				}
//	}
//
//	bool SetSize(GLint i, GLint j, GLint k, GLdouble size) {
//		if (i >= 0 && i < particle_Grid_Len
//			&& j >= 0 && j < particle_Grid_Len
//			&& k >= 0 && k < particle_Grid_Len) {
//			particle[i][j][k].size = size;
//			return true;
//		}
//		return false;
//	}
//
//	int GetSize(GLint i, GLint j, GLint k, GLdouble &size) {
//		if (i >= 0 && i < particle_Grid_Len
//			&& j >= 0 && j < particle_Grid_Len
//			&& k >= 0 && k < particle_Grid_Len) {
//			size = particle[i][j][k].size;
//			return true;
//		}
//		return false;
//	}
//
//	/* 设置和获取温度信息 */
//	void SetAllTemperature(GLdouble temp) {
//		for (int i = 0; i < particle_Grid_Len; i++)
//			for (int j = 0; j < particle_Grid_Len; j++)
//				for (int k = 0; k < particle_Grid_Len; k++) {
//					particle[i][j][k].temperature = temp;
//				}
//	}
//
//	bool SetTemperature(GLint i, GLint j, GLint k, GLdouble temp) {
//		if (i >= 0 && i < particle_Grid_Len
//			&& j >= 0 && j < particle_Grid_Len
//			&& k >= 0 && k < particle_Grid_Len) {
//			particle[i][j][k].temperature = temp;
//			return true;
//		}
//		return false;
//	}
//
//	int GetTemperature(GLint i, GLint j, GLint k, GLdouble &temp) {
//		if (i >= 0 && i < particle_Grid_Len
//			&& j >= 0 && j < particle_Grid_Len
//			&& k >= 0 && k < particle_Grid_Len) {
//			temp = particle[i][j][k].temperature;
//			return true;
//		}
//		return false;
//	}
//
//	void CalAvgTemperature() {
//		for (int i = 1; i < particle_Grid_Len - 1; i++)
//			for (int j = 1; j < particle_Grid_Len - 1; j++)
//				for (int k = 1; k < particle_Grid_Len - 1; k++) {
//					temp_avg += particle[i][j][k].temperature;
//				}
//		temp_avg /= (particle_Grid_Len - 2)* (particle_Grid_Len - 2)* (particle_Grid_Len - 2);
//	}
//
//	/* 设置和获取密度信息 */
//	void SetAllDensity(GLdouble density) {
//		for (int i = 0; i < particle_Grid_Len; i++)
//			for (int j = 0; j < particle_Grid_Len; j++)
//				for (int k = 0; k < particle_Grid_Len; k++) {
//					particle[i][j][k].density = density;
//				}
//	}
//
//	bool SetDensity(GLint i, GLint j, GLint k, GLdouble density) {
//		if (i >= 0 && i < particle_Grid_Len
//			&& j >= 0 && j < particle_Grid_Len
//			&& k >= 0 && k < particle_Grid_Len) {
//			particle[i][j][k].density = density;
//			return true;
//		}
//		return false;
//	}
//
//	int GetDensity(GLint i, GLint j, GLint k, GLdouble &density) {
//		if (i >= 0 && i < particle_Grid_Len
//			&& j >= 0 && j < particle_Grid_Len
//			&& k >= 0 && k < particle_Grid_Len) {
//			density = particle[i][j][k].density;
//			return true;
//		}
//		return false;
//	}
//
//	/* 获取粒子数组地址 */
//	//Particle *GetParticle() { 
//	//	return particle; 
//	//}
//
//	/* 获得粒子的数目 */
//	//int GetParticleCnt() { return particle_Cnt; }
//
//	/* 设置粒子的所有属性 */
//	int SetAllData(int i, int j, int k,                /* 下标 */
//		GLint r, GLint g, GLint b,                /* 颜色 */
//		GLdouble vx, GLdouble vy, GLdouble vz,    /* 速度 */
//		GLdouble size,                              /* 大小 */
//		GLdouble duration,                          /* 持续时间 */
//		GLdouble temperature,                       /* 温度 */
//		GLdouble density                            /* 密度 */
//		) {
//		if (i >= 0 && i < particle_Grid_Len
//			&& j >= 0 && j < particle_Grid_Len
//			&& k >= 0 && k < particle_Grid_Len) {
//			particle[i][j][k].r = r;
//			particle[i][j][k].g = g;
//			particle[i][j][k].b = b;
//			particle[i][j][k].vx = vx;
//			particle[i][j][k].vy = vy;
//			particle[i][j][k].vz = vz;
//			particle[i][j][k].size = size;
//			particle[i][j][k].duration = duration;
//			particle[i][j][k].temperature = temperature;
//			particle[i][j][k].density = density;
//			return true;
//		}
//		return false;
//	}
//
//	/* 获得粒子所有的属性 */
//	int GetAllData(int i, int j, int k,              /* 下标 */
//		GLint &r, GLint &g, GLint &b,                /* 颜色 */
//		GLdouble &vx, GLdouble &vy, GLdouble &vz,    /* 速度 */
//		GLdouble &size,                              /* 大小 */
//		GLdouble &duration,                          /* 持续时间 */
//		GLdouble &temperature,                       /* 温度 */
//		GLdouble &density                            /* 密度 */
//		) {
//
//		if (i >= 0 && i < particle_Grid_Len
//			&& j >= 0 && j < particle_Grid_Len
//			&& k >= 0 && k < particle_Grid_Len) {
//			r = particle[i][j][k].r;
//			g = particle[i][j][k].g;
//			b = particle[i][j][k].b;
//			vx = particle[i][j][k].vx;
//			vy = particle[i][j][k].vy;
//			vz = particle[i][j][k].vz;
//			size = particle[i][j][k].size;
//			duration = particle[i][j][k].duration;
//			temperature = particle[i][j][k].temperature;
//			density = particle[i][j][k].density;
//			return true;
//		}
//		return false;
//	}
//
//	void initAllData() {
//		for (int i = 0; i < particle_Grid_Len; i++)
//			for (int j = 0; j < particle_Grid_Len; j++)
//				for (int k = 0; k < particle_Grid_Len; k++) {
//					SetAllData(i, j, k, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
//					for (int m = 0; m < 3; m++)
//						particle[i][j][k].force[m] = 0;
//				}
//	}
//
//	bool calculateVorticity(int i, int j, int k, double vort[3]) {
//		if (i > 0 && i < particle_Grid_Len - 1
//			&& j > 0 && j < particle_Grid_Len - 1
//			&& k > 0 && k < particle_Grid_Len - 1) {
//			vort[0] = (particle[i][j + 1][k - 1].vz + particle[i][j + 1][k].vz
//				- particle[i][j - 1][k - 1].vz - particle[i][j - 1][k].vz
//				- particle[i][j - 1][k + 1].vy - particle[i][j][k + 1].vy
//				+ particle[i][j - 1][k - 1].vy + particle[i][j][k - 1].vy) / (4 * h);
//
//			vort[1] = (particle[i - 1][j][k + 1].vx + particle[i][j][k + 1].vx
//				- particle[i - 1][j][k - 1].vx - particle[i][j][k - 1].vx
//				- particle[i + 1][j][k - 1].vz - particle[i + 1][j][k].vz
//				+ particle[i - 1][j][k - 1].vz + particle[i - 1][j][k].vz) / (4 * h);
//
//			vort[2] = (particle[i + 1][j][k - 1].vy + particle[i + 1][j][k].vy
//				- particle[i - 1][j][k - 1].vy - particle[i - 1][j][k].vy
//				- particle[i - 1][j + 1][k].vx - particle[i][j + 1][k].vx
//				+ particle[i - 1][j - 1][k].vx + particle[i][j - 1][k].vx) / (4 * h);
//			return true;
//		}
//		return false;
//	}
//
//	bool calculateForce(int i,int j,int k) {
//		if (i > 0 && i < particle_Grid_Len - 1
//			&& j > 0 && j < particle_Grid_Len - 1
//			&& k > 0 && k < particle_Grid_Len - 1) {
//			double alpha = 0.0001;
//			double beta = 0.0001;
//			double eposilo = 0.0001;
//			double array_z[3] = { 0,0,1 };
//
//			double x_next_vort[3] = {};
//			double y_next_vort[3] = {};
//			double z_next_vort[3] = {};
//			double vorticity[3] = {};
//			double nita[3] = {};
//			double vort_loca[3] = {};
//			//double force[3] = {};
//
//
//			for (int m = 0; m < 3; m++) {
//				particle[i][j][k].force[m] = 0;
//				particle[i][j][k].force[m] += beta * (particle[i][j][k].temperature - temp_avg) * array_z[m];
//				particle[i][j][k].force[m] -= alpha * particle[i][j][k].density * array_z[m];
//			}
//
//			calculateVorticity(i, j, k, vorticity);
//			calculateVorticity(i + 1, j, k, x_next_vort);
//			calculateVorticity(i, j + 1, k, y_next_vort);
//			calculateVorticity(i, j, k + 1, z_next_vort);
//			
//			double omega = sqrt(vorticity[0] * vorticity[0] +
//				vorticity[1] * vorticity[1] + vorticity[2] * vorticity[2]);
//
//			nita[0] = (sqrt(x_next_vort[0] * x_next_vort[0] + x_next_vort[1] * x_next_vort[1]
//				+x_next_vort[2] * x_next_vort[2]) - omega) / h;
//			nita[1] = (sqrt(y_next_vort[0] * y_next_vort[0] + y_next_vort[1] * y_next_vort[1]
//				+ y_next_vort[2] * y_next_vort[2]) - omega) / h;
//			nita[2] = (sqrt(z_next_vort[0] * z_next_vort[0] + z_next_vort[1] * z_next_vort[1]
//				+ z_next_vort[2] * z_next_vort[2]) - omega) / h;
//
//			double nita_len = sqrt(nita[0] * nita[0] + nita[1] * nita[1] + nita[2] * nita[2]);
//			
//			if (nita_len != 0) {
//				for (int m = 0; m < 3; m++)
//					vort_loca[m] = nita[m] / nita_len;
//			}
//
//			particle[i][j][k].force[0] += eposilo * h * (vort_loca[1] * vorticity[2] - vort_loca[2] * vorticity[1]);
//			particle[i][j][k].force[1] += eposilo * h * (vort_loca[2] * vorticity[0] - vort_loca[0] * vorticity[2]);
//			particle[i][j][k].force[2] += eposilo * h * (vort_loca[0] * vorticity[1] - vort_loca[1] * vorticity[0]);
//
//			return true;
//		}
//		return false;
//	}
//
//}minus_smoke,zero_smoke,plus_smoke;


double vx, vy, vz, vStarX,vStarY,vStarZ,size, duration, temperature, density;
int r, g, b;
double deltaTime = 0.01;
Vec3 advectionTerm;

void InitSmoke() {

	/* 初始化颜色 */
	r = 192;
	g = 192;
	b = 192;
	old_smoke.SetAllColor(r, g, b);

	old_smoke.SetAllVelocity(0,0,0);
	old_smoke.SetAllDuration(100);
	old_smoke.SetAllSize(0.08f);
	old_smoke.SetAllTemperature(0);
	old_smoke.CalAvgTemperature();

	/* 初始化密度 */
	for (int i = 0; i < particle_Grid_Len; i++)
		for (int j = 0; j < particle_Grid_Len; j++)
			for (int k = 0; k < particle_Grid_Len; k++) {
				double x = 4 * ((double)i / particle_Grid_Len) - 2;
				double y = 4 * ((double)j / particle_Grid_Len) - 2;
				double z = 4 * ((double)k / particle_Grid_Len) - 2;
				if (x * x + y * y + z * z < 1)
					old_smoke.SetDensity(i, j, k, 1 - x*x - y*y - z*z);
				else 
					old_smoke.SetDensity(i, j, k, 0);
			}
}

static const float vertex_list[][3] =
{
	-0.5f, -0.5f, -0.5f,
	0.5f, -0.5f, -0.5f,
	-0.5f, 0.5f, -0.5f,
	0.5f, 0.5f, -0.5f,
	-0.5f, -0.5f, 0.5f,
	0.5f, -0.5f, 0.5f,
	-0.5f, 0.5f, 0.5f,
	0.5f, 0.5f, 0.5f,
};

// 将要使用的顶点的序号保存到一个数组里面   

static const GLint index_list[][2] =
{
	{ 0, 1 },
	{ 2, 3 },
	{ 4, 5 },
	{ 6, 7 },
	{ 0, 2 },
	{ 1, 3 },
	{ 4, 6 },
	{ 5, 7 },
	{ 0, 4 },
	{ 1, 5 },
	{ 7, 3 },
	{ 2, 6 }
};

// 绘制立方体  

void DrawEdge(void)
{
	int i, j;
	glColor3b(0, 0, 0);
	glBegin(GL_LINES);
	for (i = 0; i<12; ++i) // 12 条线段  

	{
		for (j = 0; j<2; ++j) // 每条线段 2个顶点  

		{
			glVertex3fv(vertex_list[index_list[i][j]]);
		}
	}
	glEnd();
}

//static int count = 0;

void DrawSmoke() {

		//count++;
		old_smoke.CalAvgTemperature();

		for (int i = 1; i < particle_Grid_Len - 1; i++)
			for (int j = 1; j < particle_Grid_Len - 1; j++)
				for (int k = 1; k < particle_Grid_Len - 1; k++) {

					old_smoke.GetAllData(i, j, k, r, g, b,
						vx, vy, vz, vStarX,vStarY,vStarZ,size, duration, temperature, density,advectionTerm);

					glColor4ub(r, g, b, 255);
					//glNormal3f(0.0f, 0.0f, 1.0f);

					//double x_location = 2 * ((double)i / particle_Grid_Len) - 1;
					//double y_location = 2 * ((double)j / particle_Grid_Len) - 1;
					//double z_location = 2 * ((double)k / particle_Grid_Len) - 1;

					double x_location = ((double)i / particle_Grid_Len) - 0.5;
					double y_location = ((double)j / particle_Grid_Len) - 0.5;
					double z_location = ((double)k / particle_Grid_Len) - 0.5;

					//x_location = x_location/count;
					//y_location = y_location/count;
					//z_location = z_location/count;

					//double x_location = ((double)i / particle_Grid_Len);
					//double y_location = ((double)j / particle_Grid_Len);
					//double z_location = ((double)k / particle_Grid_Len);

					//glTranslated(x_location, y_location, z_location);
					//GLUquadricObj* mySphere = gluNewQuadric();
					//gluSphere(mySphere, size, 32, 16);

					if (density > 0) {
						glBegin(GL_QUADS);
						glVertex3d(x_location - density * size, y_location - density * size, z_location);
						glVertex3d(x_location - density * size, y_location + density * size, z_location);
						glVertex3d(x_location + density * size, y_location + density * size, z_location);
						glVertex3d(x_location + density * size, y_location - density * size, z_location);
						glEnd();
						//glBegin(GL_QUADS);
						//glVertex3d(x_location - size, y_location - size, z_location);
						//glVertex3d(x_location - size, y_location + size, z_location);
						//glVertex3d(x_location + size, y_location + size, z_location);
						//glVertex3d(x_location + size, y_location - size, z_location);
						//glEnd();
						//glFlush();
					}

					old_smoke.calculateForce(i, j, k);
				}
				
		
		glFlush();
		for (int i = 1; i < particle_Grid_Len - 1; i++)
			for (int j = 1; j < particle_Grid_Len - 1; j++)
				for (int k = 1; k < particle_Grid_Len - 1; k++) {
					old_smoke.GetAllData(i, j, k, r, g, b,
						vx, vy, vz, vStarX,vStarY,vStarZ,size, duration, temperature, density, advectionTerm);
					vx += deltaTime * (old_smoke.particle[i][j][k].force[0]
						+ old_smoke.particle[i + 1][j][k].force[0]) / 2.0;
					vy += deltaTime * (old_smoke.particle[i][j][k].force[1]
						+ old_smoke.particle[i][j + 1][k].force[1]) / 2.0;
					vz += deltaTime * (old_smoke.particle[i][j][k].force[2]
						+ old_smoke.particle[i][j][k + 1].force[2]) / 2.0;
					old_smoke.SetAllData(i, j, k, r, g, b,
						vx, vy, vz, vStarX, vStarY, vStarZ, size, duration, temperature, density,advectionTerm);
				}

		semiLagrangeCalc(old_smoke, gen_smoke);

		for (int i = 1; i < particle_Grid_Len - 1; i++)
			for (int j = 1; j < particle_Grid_Len - 1; j++)
				for (int k = 1; k < particle_Grid_Len - 1; k++) {
					gen_smoke.GetAllData(i, j, k, r, g, b,
						vx, vy, vz, vStarX, vStarY, vStarZ, size, duration, temperature, density,advectionTerm);
					old_smoke.SetAllData(i, j, k, r, g, b,
						vx, vy, vz, vStarX, vStarY, vStarZ, size, duration, temperature, density,advectionTerm);
				}
		gen_smoke.initAllData();
		glutPostRedisplay();
}

void myInit() {

	//glMatrixMode(GL_PROJECTION);
	//glLoadIdentity();
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
	glClearDepth(1.0f);                   //指定深度缓冲区的清除值(即将深度缓冲区里的值设置为1)
	//glDepthFunc(GL_LEQUAL);               //指定用于深度缓冲比较值(即新进像素深度值与原来的1比较，<=则通过，否则丢弃)
	glEnable(GL_DEPTH_TEST);	
	glDepthFunc(GL_ALWAYS);
	glShadeModel(GL_SMOOTH);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

	old_smoke.CreateGrid(particle_Grid_Len);
	old_smoke.initAllData();
	gen_smoke.CreateGrid(particle_Grid_Len);
	gen_smoke.initAllData();
	
	InitSmoke();
	//glFlush();

}

void display() {

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//glLoadIdentity();
	DrawEdge();
	DrawSmoke();

	//DrawEdge();

	glFlush();
}

void ChangeSize(int width, int height)
{
	GLfloat aspect = (GLfloat)width / (GLfloat)height;
	GLfloat nRange = 100.0f;
	/* 设置视口大小为整个窗口大小 */
	glViewport(0, 0, width, height);
	/* 单位化投影矩阵 */
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//glFrustum(-0.5F, 0.5F, -0.5F, 0.5F, 1.0F, 3.0F);
	/* 设置正确的投影矩阵 */
	//gluPerspective(45.0f, (GLfloat)width / (GLfloat)height, 0.1f, 1000.0f);
	//设置三维投影区  
	//if (width <= height)
	//{
	//	glOrtho(-nRange, nRange, -nRange * aspect, nRange * aspect, -nRange, nRange);
	//}
	//else
	//{
	//	glOrtho(-nRange, nRange, -nRange / aspect, nRange / aspect, -nRange, nRange);
	//}
	/* 设置模型视图矩阵 */
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glRotatef(20.0F, 1.0F, 0.0F, 0.0F);
	glRotatef(-30.0F, 0.0F, 1.0F, 0.0F);
	//gluLookAt(0.0f, 0.0f, 10.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
}


int main(int argc, char** argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_RGB);
	glutInitWindowSize(800, 800);
	//glutInitWindowPosition(0, 0);
	glutCreateWindow("Smoke Simulation");
	myInit();
	glutReshapeFunc(ChangeSize);
	glutDisplayFunc(display);
	glutMainLoop();
}



