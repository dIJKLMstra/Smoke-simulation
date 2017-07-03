#include <windows.h>
#include <time.h>
#include <stdlib.h>
#include <gl\GL.h>
#include <gl\GLU.h>
#include <gl\glut.h>
#include <GL\GLAUX.H>
#include <iostream>

#include "semiLag.hpp"
#include "semiLag.cpp"
#include "vector.hpp"
#include "vector.cpp"

#pragma comment(lib, "glaux.lib")

/* 使用到的一些参数 */

/* 用来保存当前粒子属性的一些变量 */
double vx, vy, vz, vStarX, vStarY, vStarZ;
double size, duration, temperature, density;
int r, g, b;

double extinction = 1;					/* 消光系数 */
double e = 2.718281828459;
double deltaTime = 0.01;				/* 每一帧的时间增量 */
Vec3 advectionTerm;
GLuint	texture[1];

/* 加载纹理 */
AUX_RGBImageRec *LoadBMP(char *Filename)				
{
	FILE *File = NULL;								
	if (!Filename)									
	{
		return NULL;							
	}
	File = fopen(Filename, "r");						
	if (File)										
	{
		fclose(File);								
		return auxDIBImageLoad(Filename);			
	}
	return NULL;									
}

/* 将图片转成纹理 */
int LoadGLTextures()								
{
	int Status = FALSE;							
	AUX_RGBImageRec *TextureImage[1];				
	memset(TextureImage, 0, sizeof(void *) * 1);	

	if (TextureImage[0] = LoadBMP("smoke.bmp"))	    
	{
		Status = TRUE;								
		glGenTextures(1, &texture[0]);				

		glBindTexture(GL_TEXTURE_2D, texture[0]);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexImage2D(GL_TEXTURE_2D, 0, 3, TextureImage[0]->sizeX, 
			TextureImage[0]->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, TextureImage[0]->data);
	}

	if (TextureImage[0])						
	{
		if (TextureImage[0]->data)					
		{
			free(TextureImage[0]->data);			
		}
		free(TextureImage[0]);					
	}
	return Status;									
}

/* 对烟雾的初始化 */
void InitSmoke() {

	srand(time(NULL));

	/* 初始化颜色 */
	r = 100;
	g = 100;
	b = 100;
	//r = rand() % 255;
	//g = rand() % 255;
	//b = rand() % 255;
	old_smoke.SetAllColor(r, g, b);

	/* 初始化持续时间 尺寸 温度 */
	//old_smoke.SetAllVelocity(0,0,0);
	old_smoke.SetAllDuration(100);
	old_smoke.SetAllSize(1.5 / particle_Grid_Len);
	old_smoke.SetAllTemperature(0);
	//old_smoke.CalAvgTemperature();

	/* 初始化密度与速度 */
	for (int i = 0; i < particle_Grid_Len; i++)
		for (int j = 0; j < particle_Grid_Len; j++)
			for (int k = 0; k < particle_Grid_Len; k++) {
				//double x = 4 * ((double)i / particle_Grid_Len) - 2;
				//double y = 4 * ((double)j / particle_Grid_Len) - 2;
				//double z = 4 * ((double)k / particle_Grid_Len) - 2;
				//if (x * x + y * y + z * z < 1)
				//	old_smoke.SetDensity(i, j, k, 1 - x*x - y*y - z*z);
				//else 
				//	old_smoke.SetDensity(i, j, k, 0);
				double x = 20 * ((double)i / particle_Grid_Len) - 10;
				double y = 5 * ((double)j / particle_Grid_Len);
				double z = 20 * ((double)k / particle_Grid_Len) - 10;

				/* 在下表面中心位置产生烟雾 */
				if (x * x + z * z < 1 && y < 1) {
					old_smoke.SetDensity(i, j, k, 
						(1 - x * x - z * z) * (1 - y) * (rand() % 1000 / 1000.0));
					old_smoke.SetVelocity(i, j, k, (rand() % 1000 / 500.0 - 1), 
						(rand() % 1000 / 500.0 - 1), (rand() % 1000 / 500.0 - 1));
					//old_smoke.SetTemperature(i, j, k, (rand() % 1000 / 100000.0));
					//old_smoke.SetVelocity(i, j, k, 0,0,0);
				}
				else {
					old_smoke.SetDensity(i, j, k, 0);
					old_smoke.SetVelocity(i, j, k, 0,0,0);
					//old_smoke.SetTemperature(i, j, k, 0);
				}
				/* 初始化V*与V一致 */
				old_smoke.SetVstar(i,j,k);
			}

	old_smoke.CalAvgTemperature();

}

/* 外边框的顶点 */
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
  
/* 外边框顶点序号 */
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

/* 绘制外边界 */
void DrawEdge(void)
{

	/* 边框材质 */
	GLfloat mat_ambient[] = { 0.0f,0.0f,0.0f,1.0f };

	/* 漫反射性质 */
	GLfloat mat_diffuse[] = { 1.0f,1.0f,1.0f,1.0f };

	/* 镜面反射性质 */
	GLfloat mat_specular[] = { 1.0f,1.0f,1.0f,1.0f };

	/* 镜面反射的光亮度 */
	GLfloat mat_shininess = 50.0;

	/* 将以上材质定义应用 */
	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialf(GL_FRONT, GL_SHININESS, mat_shininess);

	int i, j;
	//glColor3b(0, 0, 0);
	glBegin(GL_LINES);

	/* 12条线段 */
	for (i = 0; i<12; ++i) 
	{
		for (j = 0; j<2; ++j)  
		{
			glVertex3fv(vertex_list[index_list[i][j]]);
		}
	}

	glEnd();

}

/* 绘制物体边界(正方体) */
void DrawBlock(double dSize) {

	/* 正方体材质 */
	GLfloat mat_ambient[] = { 0.2f,0.3f,0.4f,1.0f };

	/* 漫反射性质 */
	GLfloat mat_diffuse[] = { 1.0f,1.0f,1.0f,1.0f };

	/* 镜面反射性质 */
	GLfloat mat_specular[] = { 1.0f,1.0f,1.0f,1.0f };

	/* 镜面反射的光亮度 */
	GLfloat mat_shininess = 30.0;

	/* 将以上材质定义应用 */
	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialf(GL_FRONT, GL_SHININESS, mat_shininess);

	double size = dSize * 0.5;
	
	/* 上表面 */
	glBegin(GL_QUADS);
	//glColor4f(1.0, 0.0, 0.0, 1.0);
	glNormal3d(0.0, 0.0, 1.0);   
	glVertex3d(dSize, dSize, dSize);
	glVertex3d(-dSize, dSize, dSize);
	glVertex3d(-dSize, -dSize, dSize);
	glVertex3d(dSize, -dSize, dSize);
	glEnd();

	/* 下表面 */
	//glColor4f(0.0, 1.0, 0.0, 1.0);
	glBegin(GL_QUADS);
	glNormal3d(0.0, 0.0, -1.0);
	glVertex3d(dSize, dSize, -dSize);
	glVertex3d(-dSize, dSize, -dSize);
	glVertex3d(-dSize, -dSize, -dSize);
	glVertex3d(dSize, -dSize, -dSize);
	glEnd();

	/* 前表面 */
	//glColor4f(1.0, 1.0, 0.0, 1.0);
	glBegin(GL_QUADS);
	glNormal3d(1.0, 0.0, 0.0); 
	glVertex3d(dSize, dSize, dSize);
	glVertex3d(dSize, -dSize, dSize);
	glVertex3d(dSize, -dSize, -dSize);
	glVertex3d(dSize, dSize, -dSize);
	glEnd();

	/* 后表面 */
	//glColor4f(0.0, 0.0, 1.0, 1.0);
	glBegin(GL_QUADS);
	glNormal3d(-1.0, 0.0, 0.0);
	glVertex3d(-dSize, dSize, dSize);
	glVertex3d(-dSize, dSize, -dSize);
	glVertex3d(-dSize, -dSize, -dSize);
	glVertex3d(-dSize, -dSize, dSize);
	glEnd();

	/* 左表面 */  
	//glColor4f(0.0, 1.0, 0.5, 1.0);
	glBegin(GL_QUADS);
	glNormal3d(0.0, -1.0, 0.0);  
	glVertex3d(dSize, -dSize, dSize);
	glVertex3d(dSize, -dSize, -dSize);
	glVertex3d(-dSize, -dSize, -dSize);
	glVertex3d(-dSize, -dSize, dSize);
	glEnd();

	/* 右表面 */  
	//glColor4f(0.5, 1.0, 0.5, 1.0);
	glBegin(GL_QUADS);
	glNormal3d(0.0, 1.0, 0.0);  
	glVertex3d(dSize, dSize, dSize);
	glVertex3d(dSize, dSize, -dSize);
	glVertex3d(-dSize, dSize, -dSize);
	glVertex3d(-dSize, dSize, dSize);  
	glEnd();

}

/* 设置光照 */
void setLight()
{
	/* 光源位置 */
	GLfloat light_position[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	
	/* 多次反射遗留在环境中的光线强度 */
	GLfloat light_ambient[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	
	/* 漫反射后光线强度 */
	GLfloat light_diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	
	/* 镜面反射后光照强度 */
	GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };

	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);

	/* 在后面的渲染中使用0号光照 */
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

}

/* 绘制烟雾 */
void DrawSmoke() {
		
		/* 先计算当前网格的平均温度 */
		old_smoke.CalAvgTemperature();

		for (int i = 1; i < particle_Grid_Len - 1; i++)
			for (int j = 1; j < particle_Grid_Len - 1; j++)
				for (int k = 1; k < particle_Grid_Len - 1; k++) {

					/* 如果是物体边界 不用计算外力 */
					if (old_smoke.particle[i][j][k].isOccupied)
						continue;

					/* 去除一些边缘烟雾 */
					old_smoke.isloateParticle(i,j,k);

					/* 获取现有烟雾网格的所有数据 */
					old_smoke.GetAllData(i, j, k, r, g, b,vx, vy, vz, 
						vStarX, vStarY, vStarZ, size, duration, temperature, density, advectionTerm);

					/* 根据当前粒子的密度计算该粒子的透明度 */
					double t_vox = 0;
					if (density > 0)
						t_vox = pow(e, -(extinction * h / density));
					double color = 0.8 * (1 - t_vox);
					//GLfloat mat_ambient[] = { r / 255.0, g / 255.0, b / 255.0, density};
					//GLfloat mat_diffuse[] = { r / 255.0, g / 255.0, b / 255.0, density};
					//GLfloat mat_specular[] = { r / 255.0, g / 255.0, b / 255.0, density};
					//GLfloat mat_emission[] = { 0.0f, 0.0f, 0.0f, density};
					GLfloat mat_ambient[] = { r / 255.0, g / 255.0, b / 255.0, t_vox};
					GLfloat mat_diffuse[] = { r / 255.0, g / 255.0, b / 255.0, 0.2};
					GLfloat mat_specular[] = { r / 255.0, g / 255.0, b / 255.0,t_vox};
					GLfloat mat_emission[] = { 0.0f, 0.0f, 0.0f, t_vox};

					GLfloat mat_shininess = 20.0f;

					glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
					glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
					glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
					glMaterialfv(GL_FRONT, GL_EMISSION, mat_emission);
					glMaterialf(GL_FRONT, GL_SHININESS, mat_shininess);

					glColor4ub(r, g, b, density * 255);
					//glNormal3f(0.0f, 0.0f, 1.0f);

					//double x_location = 2 * ((double)i / particle_Grid_Len) - 1;
					//double y_location = 2 * ((double)j / particle_Grid_Len) - 1;
					//double z_location = 2 * ((double)k / particle_Grid_Len) - 1;

					double x_location = ((double)i / particle_Grid_Len) - 0.5;
					double y_location = ((double)j / particle_Grid_Len) - 0.5;
					double z_location = ((double)k / particle_Grid_Len) - 0.5;

					//double x_location = ((double)i / particle_Grid_Len);
					//double y_location = ((double)j / particle_Grid_Len);
					//double z_location = ((double)k / particle_Grid_Len);

					/* 
					 * 绘制烟雾
					 * 绘制条件为该粒子密度大于0
					 * 持续时间也大于0
					 */
					if (density > 0 && duration > 0 
						&& !old_smoke.particle[i][j][k].border) {

						/* 调整绘图原点 */
						glTranslated(x_location, y_location, z_location);

						/* 每个烟雾粒子以球体的形式绘制 */
						glutSolidSphere(size, 10, 10);

						/* 恢复绘图原点 */
						glTranslated(-x_location, -y_location, -z_location);
						
					}

					/*
					 * 曾经使用的纹理贴图来绘制烟雾
					 * 由于效果不好并没有使用 
					 */
					/*if (density > 0) {
						glBegin(GL_QUADS);
						glTexCoord2f(0.0f, 0.0f); 
						glVertex3d(x_location - density * size, y_location - density * size, z_location);
						glTexCoord2f(1.0f, 0.0f);
						glVertex3d(x_location - density * size, y_location + density * size, z_location);
						glTexCoord2f(1.0f, 1.0f);
						glVertex3d(x_location + density * size, y_location + density * size, z_location);
						glTexCoord2f(0.0f, 1.0f);
						glVertex3d(x_location + density * size, y_location - density * size, z_location);
						glEnd();
						glFlush();
					}*/

					/* 计算这个粒子所受外力 */
					old_smoke.calculateForce(i, j, k);
				}
				
		glFlush();

		for (int i = 1; i < particle_Grid_Len - 1; i++)
			for (int j = 1; j < particle_Grid_Len - 1; j++)
				for (int k = 1; k < particle_Grid_Len - 1; k++) {

					old_smoke.GetAllData(i, j, k, r, g, b,vx, vy, vz,
						 vStarX,vStarY,vStarZ,size, duration, temperature, density, advectionTerm);

					/* 根据所受外力更新速度场 */
					vx += deltaTime * (old_smoke.particle[i][j][k].force[0]
						+ old_smoke.particle[i + 1][j][k].force[0]) / 2.0;
					vy += deltaTime * (old_smoke.particle[i][j][k].force[1]
						+ old_smoke.particle[i][j + 1][k].force[1]) / 2.0;
					vz += deltaTime * (old_smoke.particle[i][j][k].force[2]
						+ old_smoke.particle[i][j][k + 1].force[2]) / 2.0;
					//vx += 1;
					//if ((vx > 0 && i == 1) ||(vx < 0 && i == particle_Grid_Len - 2))
					//	vx = 0;
					//if ((vy > 0 && j == 1) || (vy < 0 && j == particle_Grid_Len - 2))
					//	vy = 0;
					//if ((vz > 0 && k == 1) || (vz < 0 && k == particle_Grid_Len - 2))
					//	vz = 0;

					/* 边界处理 */
					if (vx > 0 && old_smoke.particle[i - 1][j][k].border ||
						vx < 0 && old_smoke.particle[i + 1][j][k].border)
						vx = 0;
					if (vy > 0 && old_smoke.particle[i][j - 1][k].border ||
						vy < 0 && old_smoke.particle[i][j + 1][k].border)
						vy = 0;
					if (vz > 0 && old_smoke.particle[i][j][k - 1].border ||
						vz < 0 && old_smoke.particle[i][j][k + 1].border)
						vz = 0;

					
					//old_smoke.SetAllData(i, j, k, r, g, b,
					//	vx, vy, vz, vStarX, vStarY, vStarZ, size, duration, temperature, density,advectionTerm);
					old_smoke.SetVelocity(i, j, k, vx, vy, vz);
					gen_smoke.SetDuration(i, j, k, duration);
				}

		/* 
		 * 用半拉格朗日方法与old_smoke的信息
		 * 更新当前烟雾网格的速度、温度、密度等属性
		 */
		semiLagrangeCalc(old_smoke, gen_smoke);

		/* 烟雾生成区产生新的烟雾 */
		gen_smoke.replenishDensity();

		/* 将gen_smoke中的数据传递给old_smoke中 */
		for (int i = 1; i < particle_Grid_Len - 1; i++)
			for (int j = 1; j < particle_Grid_Len - 1; j++)
				for (int k = 1; k < particle_Grid_Len - 1; k++) {
					gen_smoke.GetAllData(i, j, k, r, g, b,vx, vy, vz, 
						vStarX, vStarY, vStarZ, size, duration, temperature, density,advectionTerm);
					//if ((vx > 0 && i == 1) || (vx < 0 && i == particle_Grid_Len - 2))
					//	vx = 0;
					//if ((vy > 0 && j == 1) || (vy < 0 && j == particle_Grid_Len - 2))
					//	vy = 0;
					//if ((vz > 0 && k == 1) || (vz < 0 && k == particle_Grid_Len - 2))
					//	vz = 0;

					/* 边界处理 */
					if (vx > 0 && gen_smoke.particle[i - 1][j][k].border ||
						vx < 0 && gen_smoke.particle[i + 1][j][k].border)
						vx = 0;
					if (vy > 0 && gen_smoke.particle[i][j - 1][k].border ||
						vy < 0 && gen_smoke.particle[i][j + 1][k].border)
						vy = 0;
					if (vz > 0 && gen_smoke.particle[i][j][k - 1].border ||
						vz < 0 && gen_smoke.particle[i][j][k + 1].border)
						vz = 0;

					/* 设置持续时间与密度相关 */
					duration = 800 * density - 1;

					old_smoke.SetAllData(i, j, k, r, g, b,vx, vy, vz,
						 vStarX, vStarY, vStarZ, size, duration, temperature, density,advectionTerm);
				}

		gen_smoke.initAllData();

		/* 再次调用display函数 */
		glutPostRedisplay();

}

/* 绘图相关初始化 */
void myInit() {

	/* 设置背景颜色为白色 */
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
	       
	/* 开启深度测试 */
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL); 	
	glClearDepth(1.0f);
	//glDepthFunc(GL_ALWAYS);
	//glDepthMask(GL_FALSE);
	
	/* 设置光滑着色模式 */
	glShadeModel(GL_SMOOTH);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

	/* 启用二维纹理 */
	glEnable(GL_TEXTURE_2D);							
	glBindTexture(GL_TEXTURE_2D, texture[0]);			

	/* 启用色彩混合 */
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	/* 设置光照 */
	setLight();
	
	/* 两个烟雾对象的初始化 */
	old_smoke.CreateGrid(particle_Grid_Len);
	old_smoke.initAllData();
	gen_smoke.CreateGrid(particle_Grid_Len);
	gen_smoke.initAllData();
	
	InitSmoke();

}

void display() {

	/* 清空颜色缓存与深度缓存 */
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	DrawBlock(0.15);

	DrawEdge();

	DrawSmoke();

	glFlush();
}

/* 窗口大小改变时调用 */
void ChangeSize(int width, int height)
{
	GLfloat aspect = (GLfloat)width / (GLfloat)height;
	GLfloat nRange = 100.0f;

	/* 设置视口大小为整个窗口大小 */
	glViewport(0, 0, width, height);

	/* 单位化投影矩阵 */
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	
	/* 设置投影矩阵 */
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

	/* 对显示图像进行适当旋转 */
	//glRotatef(20.0F, 1.0F, 0.0F, 0.0F);
	//glRotatef(-30.0F, 0.0F, 1.0F, 0.0F);
	glRotatef(-20.0F, 1.0F, 0.0F, 0.0F);
	glRotatef(20.0F, 0.0F, 1.0F, 0.0F);
	//glRotatef(90.0F, 1.0F, 0.0F, 0.0F);
	//glRotatef(5.0F, 0.0F, 1.0F, 0.0F);

	//gluLookAt(0.0f, 0.5f, 0.5f, 0.5f, 0.0f, 0.0f, 0.0f, 0.0f, -0.5f);
}

/* 处理按键消息 */
void keyboard(unsigned char key, int x, int y)
{
	switch (key) {
	
	/* 按R键重新绘制烟雾 */
	case 'r':
		InitSmoke();
		display();
		break;

	/* 按X键在X轴方向添加外力 */
	case 'x':
		for (int i = 1; i < particle_Grid_Len - 1; i++)
			for (int j = 1; j < particle_Grid_Len - 1; j++)
				for (int k = 1; k < particle_Grid_Len - 1; k++) {
					if (old_smoke.particle[i][j][k].density)
						old_smoke.particle[i][j][k].vx += 0.3;
					if (gen_smoke.particle[i][j][k].density)
						gen_smoke.particle[i][j][k].vx += 0.3;
				}
		break;

	/* 按Y键在Y轴方向添加外力 */
	case 'y':
		for (int i = 1; i < particle_Grid_Len - 1; i++)
			for (int j = 1; j < particle_Grid_Len - 1; j++)
				for (int k = 1; k < particle_Grid_Len - 1; k++) {
					if (old_smoke.particle[i][j][k].density)
						old_smoke.particle[i][j][k].vy -= 0.3;
					if (gen_smoke.particle[i][j][k].density)
						old_smoke.particle[i][j][k].vy -= 0.3;
				}
		break;

	/* 按Z键在Z轴方向添加外力 */
	case 'z':
		for (int i = 1; i < particle_Grid_Len - 1; i++)
			for (int j = 1; j < particle_Grid_Len - 1; j++)
				for (int k = 1; k < particle_Grid_Len - 1; k++) {
					if (old_smoke.particle[i][j][k].density)
						old_smoke.particle[i][j][k].vz += 0.3;
					if (gen_smoke.particle[i][j][k].density)
						old_smoke.particle[i][j][k].vz += 0.3;
				}
		break;

	/* 按ESC键为退出 */
	case 27:
		exit(0);
		break;
	}
}


int main(int argc, char** argv) {

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_RGB);

	/* 设置窗口 */
	glutInitWindowSize(800, 800);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Smoke Simulation");

	myInit();

	/* 一些回调函数 */
	glutReshapeFunc(ChangeSize);

	glutDisplayFunc(display);

	glutKeyboardFunc(keyboard);

	//glutCreateMenu(mymenu); 
	//glutAddMenuEntry("Clear Screen", 1);  
	//glutAddMenuEntry("Exit", 2);
	//glutAttachMenu(GLUT_RIGHT_BUTTON);
	
	glutMainLoop();
}



