#include <windows.h>
#include <time.h>
#include <stdlib.h>
#include <gl\GL.h>
#include <gl\GLU.h>
#include <gl\glut.h>
#include <iostream>

#include "vector.hpp"
#define particle_Grid_Len 40

const double h = 1.0 / particle_Grid_Len;

/* 粒子结构 */

struct Particle
{
	unsigned int  r, g, b;      /* 粒子的颜色 */
	double vx, vy, vz;          /* 粒子的当前速度 */
	double vStarX, vStarY, vStarZ;
	double duration;			/* 粒子持续时间 */
	double size;                /* 粒子尺寸 */
	double temperature;         /* 粒子温度 */
	double density;             /* 粒子密度 */
	bool border;                 /* 是否为边界 */
	Vec3 advectionTerm;
	double force[3];
};

/* 粒子网格类 */
class Grid
{
public:
	Particle***   particle;               /* 粒子指针 */
										  //int         particle_Cnt;         /* 粒子数目 */
	double      temp_avg;            /* 温度条件 */

									 /* 构造函数 */
	Grid() {
		particle = NULL;
		//particle_Cnt = 0;
		temp_avg = 0;
	}
	/* 析构函数 */
	~Grid() {
		for (int i = 0; i < particle_Grid_Len; i++)
			for (int j = 0; j < particle_Grid_Len; j++) {
				delete[] particle[i][j];
				particle[i][j] = NULL;
			}

		for (int i = 0; i < particle_Grid_Len; i++) {
			delete[] particle[i];
			particle[i] = NULL;
		}
		delete[] particle;
		particle = NULL;
	}

	/* 创建粒子数组 */
	void CreateGrid(long num) {
		if (particle)
			delete[] particle;

		particle = (Particle***)new Particle**[particle_Grid_Len];

		for (int i = 0; i < particle_Grid_Len; i++) {
			particle[i] = (Particle**)new Particle*[particle_Grid_Len];
		}

		for (int i = 0; i < particle_Grid_Len; i++)
			for (int j = 0; j < particle_Grid_Len; j++) {
				particle[i][j] = new Particle[particle_Grid_Len];
			}


	}

	/* 设置和获取颜色属性 */
	void SetAllColor(GLint r, GLint g, GLint b) {
		for (int i = 0; i < particle_Grid_Len; i++)
			for (int j = 0; j < particle_Grid_Len; j++)
				for (int k = 0; k < particle_Grid_Len; k++) {
					particle[i][j][k].r = r;
					particle[i][j][k].g = g;
					particle[i][j][k].b = b;
				}
	}

	bool SetColor(GLint i, GLint j, GLint k, GLint r, GLint g, GLint b) {
		if (i >= 0 && i < particle_Grid_Len
			&& j >= 0 && j < particle_Grid_Len
			&& k >= 0 && k < particle_Grid_Len) {
			particle[i][j][k].r = r;
			particle[i][j][k].g = g;
			particle[i][j][k].b = b;
			return true;
		}
		return false;
	}

	bool GetColor(GLint i, GLint j, GLint k, GLint &r, GLint &g, GLint &b) {
		if (i >= 0 && i < particle_Grid_Len
			&& j >= 0 && j < particle_Grid_Len
			&& k >= 0 && k < particle_Grid_Len) {
			r = particle[i][j][k].r;
			g = particle[i][j][k].g;
			b = particle[i][j][k].b;
			return true;
		}
		return false;
	}

	/* 设置和获取速度属性 */
	void SetAllVelocity(GLdouble vx, GLdouble vy, GLdouble vz) {
		for (int i = 0; i < particle_Grid_Len; i++)
			for (int j = 0; j < particle_Grid_Len; j++)
				for (int k = 0; k < particle_Grid_Len; k++) {
					particle[i][j][k].vx = vx;
					particle[i][j][k].vy = vy;
					particle[i][j][k].vz = vz;
				}
	}

	bool SetVelocity(GLint i, GLint j, GLint k, GLdouble vx, GLdouble vy, GLdouble vz) {
		if (i >= 0 && i < particle_Grid_Len
			&& j >= 0 && j < particle_Grid_Len
			&& k >= 0 && k < particle_Grid_Len) {
			particle[i][j][k].vx = vx;
			particle[i][j][k].vy = vy;
			particle[i][j][k].vz = vz;
			return true;
		}
		return false;
	}

	bool GetVelocity(GLint i, GLint j, GLint k, GLdouble &vx, GLdouble &vy, GLdouble &vz) {
		if (i >= 0 && i < particle_Grid_Len
			&& j >= 0 && j < particle_Grid_Len
			&& k >= 0 && k < particle_Grid_Len) {
			vx = particle[i][j][k].vx;
			vy = particle[i][j][k].vy;
			vz = particle[i][j][k].vz;
			return true;
		}
		return false;
	}

	bool SetVstar(int i, int j, int k) {
		if (i >= 0 && i < particle_Grid_Len
			&& j >= 0 && j < particle_Grid_Len
			&& k >= 0 && k < particle_Grid_Len) {
			particle[i][j][k].vStarX = particle[i][j][k].vx;
			particle[i][j][k].vStarY = particle[i][j][k].vy;
			particle[i][j][k].vStarZ = particle[i][j][k].vz;
			return true;
		}
		return false;
	}

	/* 设置和获取持续时间 */
	void SetAllDuration(GLdouble duration) {
		for (int i = 0; i < particle_Grid_Len; i++)
			for (int j = 0; j < particle_Grid_Len; j++)
				for (int k = 0; k < particle_Grid_Len; k++) {
					particle[i][j][k].duration = duration;
				}
	}

	bool SetDuration(GLint i, GLint j, GLint k, GLdouble duration) {
		if (i >= 0 && i < particle_Grid_Len
			&& j >= 0 && j < particle_Grid_Len
			&& k >= 0 && k < particle_Grid_Len) {
			particle[i][j][k].duration = duration;
			return true;
		}
		return false;
	}

	bool GetDuration(GLint i, GLint j, GLint k, GLdouble &duration) {
		if (i >= 0 && i < particle_Grid_Len
			&& j >= 0 && j < particle_Grid_Len
			&& k >= 0 && k < particle_Grid_Len) {
			duration = particle[i][j][k].duration;
			return true;
		}
		return false;
	}

	/* 设置和获取尺寸属性 */
	void SetAllSize(GLdouble size) {
		for (int i = 0; i < particle_Grid_Len; i++)
			for (int j = 0; j < particle_Grid_Len; j++)
				for (int k = 0; k < particle_Grid_Len; k++) {
					particle[i][j][k].size = size;
				}
	}

	bool SetSize(GLint i, GLint j, GLint k, GLdouble size) {
		if (i >= 0 && i < particle_Grid_Len
			&& j >= 0 && j < particle_Grid_Len
			&& k >= 0 && k < particle_Grid_Len) {
			particle[i][j][k].size = size;
			return true;
		}
		return false;
	}

	int GetSize(GLint i, GLint j, GLint k, GLdouble &size) {
		if (i >= 0 && i < particle_Grid_Len
			&& j >= 0 && j < particle_Grid_Len
			&& k >= 0 && k < particle_Grid_Len) {
			size = particle[i][j][k].size;
			return true;
		}
		return false;
	}

	/* 设置和获取温度信息 */
	void SetAllTemperature(GLdouble temp) {
		for (int i = 0; i < particle_Grid_Len; i++)
			for (int j = 0; j < particle_Grid_Len; j++)
				for (int k = 0; k < particle_Grid_Len; k++) {
					particle[i][j][k].temperature = temp;
				}
	}

	bool SetTemperature(GLint i, GLint j, GLint k, GLdouble temp) {
		if (i >= 0 && i < particle_Grid_Len
			&& j >= 0 && j < particle_Grid_Len
			&& k >= 0 && k < particle_Grid_Len) {
			particle[i][j][k].temperature = temp;
			return true;
		}
		return false;
	}

	int GetTemperature(GLint i, GLint j, GLint k, GLdouble &temp) {
		if (i >= 0 && i < particle_Grid_Len
			&& j >= 0 && j < particle_Grid_Len
			&& k >= 0 && k < particle_Grid_Len) {
			temp = particle[i][j][k].temperature;
			return true;
		}
		return false;
	}

	void CalAvgTemperature() {
		for (int i = 1; i < particle_Grid_Len - 1; i++)
			for (int j = 1; j < particle_Grid_Len - 1; j++)
				for (int k = 1; k < particle_Grid_Len - 1; k++) {
					temp_avg += particle[i][j][k].temperature;
				}
		temp_avg /= (particle_Grid_Len - 2)* (particle_Grid_Len - 2)* (particle_Grid_Len - 2);
	}

	/* 设置和获取密度信息 */
	void SetAllDensity(GLdouble density) {
		for (int i = 0; i < particle_Grid_Len; i++)
			for (int j = 0; j < particle_Grid_Len; j++)
				for (int k = 0; k < particle_Grid_Len; k++) {
					particle[i][j][k].density = density;
				}
	}

	bool SetDensity(GLint i, GLint j, GLint k, GLdouble density) {
		if (i >= 0 && i < particle_Grid_Len
			&& j >= 0 && j < particle_Grid_Len
			&& k >= 0 && k < particle_Grid_Len) {
			particle[i][j][k].density = density;
			return true;
		}
		return false;
	}

	int GetDensity(GLint i, GLint j, GLint k, GLdouble &density) {
		if (i >= 0 && i < particle_Grid_Len
			&& j >= 0 && j < particle_Grid_Len
			&& k >= 0 && k < particle_Grid_Len) {
			density = particle[i][j][k].density;
			return true;
		}
		return false;
	}

	/* 获取粒子数组地址 */
	//Particle *GetParticle() { 
	//	return particle; 
	//}

	/* 获得粒子的数目 */
	//int GetParticleCnt() { return particle_Cnt; }

	/* 设置粒子的所有属性 */
	int SetAllData(int i, int j, int k,                /* 下标 */
		GLint r, GLint g, GLint b,                /* 颜色 */
		GLdouble vx, GLdouble vy, GLdouble vz,    /* 速度 */
		GLdouble vStarX,GLdouble vStarY,GLdouble vStarZ,
		GLdouble size,                              /* 大小 */
		GLdouble duration,                          /* 持续时间 */
		GLdouble temperature,                       /* 温度 */
		GLdouble density,                            /* 密度 */
		Vec3 advectionTerm
		) {
		if (i >= 0 && i < particle_Grid_Len
			&& j >= 0 && j < particle_Grid_Len
			&& k >= 0 && k < particle_Grid_Len) {
			particle[i][j][k].r = r;
			particle[i][j][k].g = g;
			particle[i][j][k].b = b;
			particle[i][j][k].vx = vx;
			particle[i][j][k].vy = vy;
			particle[i][j][k].vz = vz;
			particle[i][j][k].vStarX = vStarX;
			particle[i][j][k].vStarY = vStarY;
			particle[i][j][k].vStarZ = vStarZ;
			particle[i][j][k].size = size;
			particle[i][j][k].duration = duration;
			particle[i][j][k].temperature = temperature;
			particle[i][j][k].density = density;
			particle[i][j][k].advectionTerm= advectionTerm;
			return true;
		}
		return false;
	}

	/* 获得粒子所有的属性 */
	int GetAllData(int i, int j, int k,              /* 下标 */
		GLint &r, GLint &g, GLint &b,                /* 颜色 */
		GLdouble &vx, GLdouble &vy, GLdouble &vz,    /* 速度 */
		GLdouble &vStarX,GLdouble &vStarY,GLdouble &vStarZ,
		GLdouble &size,                              /* 大小 */
		GLdouble &duration,                          /* 持续时间 */
		GLdouble &temperature,                       /* 温度 */
		GLdouble &density,                            /* 密度 */
		Vec3 &advectionTerm
		) {

		if (i >= 0 && i < particle_Grid_Len
			&& j >= 0 && j < particle_Grid_Len
			&& k >= 0 && k < particle_Grid_Len) {
			r = particle[i][j][k].r;
			g = particle[i][j][k].g;
			b = particle[i][j][k].b;
			vx = particle[i][j][k].vx;
			vy = particle[i][j][k].vy;
			vz = particle[i][j][k].vz;
			vStarX = particle[i][j][k].vStarX;
			vStarY = particle[i][j][k].vStarY;
			vStarZ = particle[i][j][k].vStarZ;
			size = particle[i][j][k].size;
			duration = particle[i][j][k].duration;
			temperature = particle[i][j][k].temperature;
			density = particle[i][j][k].density;
			advectionTerm = particle[i][j][k].advectionTerm;
			return true;
		}
		return false;
	}

	/* 粒子系统的初始化 */
	void initAllData() {
		Vec3 zero;
		//srand(time(NULL));
		//int r = rand() % 255;
		//int g = rand() % 255;
		//int b = rand() % 255;
		for (int i = 0; i < particle_Grid_Len; i++)
			for (int j = 0; j < particle_Grid_Len; j++)
				for (int k = 0; k < particle_Grid_Len; k++) {
					SetAllData(i, j, k,100,100,100, 0, 0, 0,
						0, 0, 0, double(1.5/particle_Grid_Len), 0, 0, 0, zero);
					
					for (int m = 0; m < 3; m++)
						particle[i][j][k].force[m] = 0;

					double x_loca = (i / particle_Grid_Len) - 0.5;
					double y_loca = (j / particle_Grid_Len) - 0.5;
					double z_loca = (k / particle_Grid_Len) - 0.5;

					/* 边界的初始化 */
					if (i == 0 || i == particle_Grid_Len - 1
						|| j == 0 || j == particle_Grid_Len - 1
						|| k == 0 || k == particle_Grid_Len - 1)
						particle[i][j][k].border = true;

					else if ((x_loca > -0.15 && i < x_loca < 0.15)
						|| (y_loca > -0.15 && y_loca < 0.15)
						|| (z_loca > -0.15 && z_loca < 0.15))
						particle[i][j][k].border = true;

					else
						particle[i][j][k].border = false;
				}
	}

	bool calculateVorticity(int i, int j, int k, double vort[3]) {
		if (i > 0 && i < particle_Grid_Len - 1
			&& j > 0 && j < particle_Grid_Len - 1
			&& k > 0 && k < particle_Grid_Len - 1) {
			vort[0] = (particle[i][j + 1][k - 1].vz + particle[i][j + 1][k].vz
				- particle[i][j - 1][k - 1].vz - particle[i][j - 1][k].vz
				- particle[i][j - 1][k + 1].vy - particle[i][j][k + 1].vy
				+ particle[i][j - 1][k - 1].vy + particle[i][j][k - 1].vy) / (4 * h);

			vort[1] = (particle[i - 1][j][k + 1].vx + particle[i][j][k + 1].vx
				- particle[i - 1][j][k - 1].vx - particle[i][j][k - 1].vx
				- particle[i + 1][j][k - 1].vz - particle[i + 1][j][k].vz
				+ particle[i - 1][j][k - 1].vz + particle[i - 1][j][k].vz) / (4 * h);

			vort[2] = (particle[i + 1][j][k - 1].vy + particle[i + 1][j][k].vy
				- particle[i - 1][j][k - 1].vy - particle[i - 1][j][k].vy
				- particle[i - 1][j + 1][k].vx - particle[i][j + 1][k].vx
				+ particle[i - 1][j - 1][k].vx + particle[i][j - 1][k].vx) / (4 * h);
			return true;
		}
		return false;
	}

	bool calculateForce(int i, int j, int k) {
		if (i > 0 && i < particle_Grid_Len - 1
			&& j > 0 && j < particle_Grid_Len - 1
			&& k > 0 && k < particle_Grid_Len - 1) {
			double alpha = particle_Grid_Len;
			double beta = particle_Grid_Len;
			double eposilo = particle_Grid_Len;
			double array_z[3] = { 0,1,0 };

			double x_next_vort[3] = {};
			double y_next_vort[3] = {};
			double z_next_vort[3] = {};
			double vorticity[3] = {};
			double nita[3] = {};
			double vort_loca[3] = {};
			//double force[3] = {};


			for (int m = 0; m < 3; m++) {
				particle[i][j][k].force[m] = 0;
				particle[i][j][k].force[m] += beta * (particle[i][j][k].temperature - temp_avg) * array_z[m];
				particle[i][j][k].force[m] -= alpha * particle[i][j][k].density * array_z[m];
			}

			calculateVorticity(i, j, k, vorticity);
			calculateVorticity(i + 1, j, k, x_next_vort);
			calculateVorticity(i, j + 1, k, y_next_vort);
			calculateVorticity(i, j, k + 1, z_next_vort);

			double omega = sqrt(vorticity[0] * vorticity[0] +
				vorticity[1] * vorticity[1] + vorticity[2] * vorticity[2]);

			nita[0] = (sqrt(x_next_vort[0] * x_next_vort[0] + x_next_vort[1] * x_next_vort[1]
				+ x_next_vort[2] * x_next_vort[2]) - omega) / h;
			nita[1] = (sqrt(y_next_vort[0] * y_next_vort[0] + y_next_vort[1] * y_next_vort[1]
				+ y_next_vort[2] * y_next_vort[2]) - omega) / h;
			nita[2] = (sqrt(z_next_vort[0] * z_next_vort[0] + z_next_vort[1] * z_next_vort[1]
				+ z_next_vort[2] * z_next_vort[2]) - omega) / h;

			double nita_len = sqrt(nita[0] * nita[0] + nita[1] * nita[1] + nita[2] * nita[2]);

			if (nita_len != 0) {
				for (int m = 0; m < 3; m++)
					vort_loca[m] = nita[m] / nita_len;
			}

			particle[i][j][k].force[0] += eposilo * h * (vort_loca[1] * vorticity[2] - vort_loca[2] * vorticity[1]);
			particle[i][j][k].force[1] += eposilo * h * (vort_loca[2] * vorticity[0] - vort_loca[0] * vorticity[2]);
			particle[i][j][k].force[2] += eposilo * h * (vort_loca[0] * vorticity[1] - vort_loca[1] * vorticity[0]);

			return true;
		}
		return false;
	}

} old_smoke, gen_smoke;
