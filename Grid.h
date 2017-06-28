#include <windows.h>
#include <gl\GL.h>
#include <gl\GLU.h>
#include <gl\glut.h>
#include <iostream>

#define particle_Grid_Len 50

const double h = 1.0 / particle_Grid_Len;

/* ���ӽṹ */

struct Particle
{
	unsigned int  r, g, b;      /* ���ӵ���ɫ */
	double vx, vy, vz;          /* ���ӵĵ�ǰ�ٶ� */
	double duration;			/* ���ӳ���ʱ�� */
	double size;                /* ���ӳߴ� */
	double temperature;         /* �����¶� */
	double density;             /* �����ܶ� */
	double force[3];
};

/* ���������� */
class Grid
{
public:
	Particle***   particle;               /* ����ָ�� */
										  //int         particle_Cnt;         /* ������Ŀ */
	double      temp_avg;            /* �¶����� */

									 /* ���캯�� */
	Grid() {
		particle = NULL;
		//particle_Cnt = 0;
		temp_avg = 0;
	}
	/* �������� */
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

	/* ������������ */
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

	/* ���úͻ�ȡ��ɫ���� */
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

	/* ���úͻ�ȡ�ٶ����� */
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

	/* ���úͻ�ȡ����ʱ�� */
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

	/* ���úͻ�ȡ�ߴ����� */
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

	/* ���úͻ�ȡ�¶���Ϣ */
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

	/* ���úͻ�ȡ�ܶ���Ϣ */
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

	/* ��ȡ���������ַ */
	//Particle *GetParticle() { 
	//	return particle; 
	//}

	/* ������ӵ���Ŀ */
	//int GetParticleCnt() { return particle_Cnt; }

	/* �������ӵ��������� */
	int SetAllData(int i, int j, int k,                /* �±� */
		GLint r, GLint g, GLint b,                /* ��ɫ */
		GLdouble vx, GLdouble vy, GLdouble vz,    /* �ٶ� */
		GLdouble size,                              /* ��С */
		GLdouble duration,                          /* ����ʱ�� */
		GLdouble temperature,                       /* �¶� */
		GLdouble density                            /* �ܶ� */
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
			particle[i][j][k].size = size;
			particle[i][j][k].duration = duration;
			particle[i][j][k].temperature = temperature;
			particle[i][j][k].density = density;
			return true;
		}
		return false;
	}

	/* ����������е����� */
	int GetAllData(int i, int j, int k,              /* �±� */
		GLint &r, GLint &g, GLint &b,                /* ��ɫ */
		GLdouble &vx, GLdouble &vy, GLdouble &vz,    /* �ٶ� */
		GLdouble &size,                              /* ��С */
		GLdouble &duration,                          /* ����ʱ�� */
		GLdouble &temperature,                       /* �¶� */
		GLdouble &density                            /* �ܶ� */
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
			size = particle[i][j][k].size;
			duration = particle[i][j][k].duration;
			temperature = particle[i][j][k].temperature;
			density = particle[i][j][k].density;
			return true;
		}
		return false;
	}

	void initAllData() {
		for (int i = 0; i < particle_Grid_Len; i++)
			for (int j = 0; j < particle_Grid_Len; j++)
				for (int k = 0; k < particle_Grid_Len; k++) {
					SetAllData(i, j, k, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
					for (int m = 0; m < 3; m++)
						particle[i][j][k].force[m] = 0;
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
			double alpha = 0.0001;
			double beta = 0.0001;
			double eposilo = 0.0001;
			double array_z[3] = { 0,0,1 };

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

}minus_smoke, zero_smoke, plus_smoke;
