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


double vx, vy, vz, vStarX,vStarY,vStarZ,size, duration, temperature, density;
int r, g, b;
double extinction = 1;
double e = 2.718281828459;
double deltaTime = 0.01;
Vec3 advectionTerm;
GLuint	texture[1];

/* 加载纹理 */
AUX_RGBImageRec *LoadBMP(char *Filename)				// Loads A Bitmap Image
{
	FILE *File = NULL;								// File Handle
	if (!Filename)									// Make Sure A Filename Was Given
	{
		return NULL;							// If Not Return NULL
	}
	File = fopen(Filename, "r");						// Check To See If The File Exists
	if (File)										// Does The File Exist?
	{
		fclose(File);								// Close The Handle
		return auxDIBImageLoad(Filename);			// Load The Bitmap And Return A Pointer
	}
	return NULL;									// If Load Failed Return NULL
}

int LoadGLTextures()									// Load Bitmap And Convert To A Texture
{
	int Status = FALSE;								// Status Indicator
	AUX_RGBImageRec *TextureImage[1];				// Create Storage Space For The Textures
	memset(TextureImage, 0, sizeof(void *) * 1);		// Set The Pointer To NULL

	if (TextureImage[0] = LoadBMP("smoke.bmp"))	    // Load Particle Texture
	{
		Status = TRUE;								// Set The Status To TRUE
		glGenTextures(1, &texture[0]);				// Create One Texture

		glBindTexture(GL_TEXTURE_2D, texture[0]);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexImage2D(GL_TEXTURE_2D, 0, 3, TextureImage[0]->sizeX, TextureImage[0]->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, TextureImage[0]->data);
	}

	if (TextureImage[0])							// If Texture Exists
	{
		if (TextureImage[0]->data)					// If Texture Image Exists
		{
			free(TextureImage[0]->data);			// Free The Texture Image Memory
		}
		free(TextureImage[0]);						// Free The Image Structure
	}
	return Status;									// Return The Status
}

/* 对烟雾的初始化 */
void InitSmoke() {

	srand(time(NULL));

	/* 初始化颜色 */
	//r = 192;
	//g = 192;
	//b = 192;
	r = rand() % 255;
	g = rand() % 255;
	b = rand() % 255;
	old_smoke.SetAllColor(r, g, b);

	//old_smoke.SetAllVelocity(0,0,0);
	old_smoke.SetAllDuration(100);
	old_smoke.SetAllSize(0.08f);
	old_smoke.SetAllTemperature(0);
	old_smoke.CalAvgTemperature();

	

	/* 初始化密度 */
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
					old_smoke.SetDensity(i, j, k, (1 - x * x - z * z) * (1 - y) * (rand() % 1000 / 1000.0));
					old_smoke.SetVelocity(i, j, k, (rand() % 1000 / 500.0 - 1), 
						(rand() % 1000 / 500.0 - 1), (rand() % 1000 / 500.0 - 1));
				}
				else {
					old_smoke.SetDensity(i, j, k, 0);
					old_smoke.SetVelocity(i, j, k, 0,0,0);
				}
				/* 初始化V*与V一致 */
				old_smoke.SetVstar(i,j,k);
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

	//定义材质
	GLfloat mat_ambient[] = { 0.0f,0.0f,0.0f,1.0f };

	//定义漫反射特性
	GLfloat mat_diffuse[] = { 1.0f,1.0f,1.0f,1.0f };

	//定义镜面反射特性
	GLfloat mat_specular[] = { 1.0f,1.0f,1.0f,1.0f };

	//定义镜面反射的光亮度
	GLfloat mat_shininess = 50.0;

	//将以上材质定义应用
	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialf(GL_FRONT, GL_SHININESS, mat_shininess);
	int i, j;
	//glColor3b(0, 0, 0);
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

void DrawBlock(double dSize) {

	//定义材质
	GLfloat mat_ambient[] = { 0.2f,0.3f,0.4f,1.0f };

	//定义漫反射特性
	GLfloat mat_diffuse[] = { 1.0f,1.0f,1.0f,1.0f };

	//定义镜面反射特性
	GLfloat mat_specular[] = { 1.0f,1.0f,1.0f,1.0f };

	//定义镜面反射的光亮度
	GLfloat mat_shininess = 30.0;

	//将以上材质定义应用
	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialf(GL_FRONT, GL_SHININESS, mat_shininess);

	double size = dSize * 0.5;
	//glEnable(GL_DEPTH_TEST);//深度测试消除隐藏面  
	glBegin(GL_QUADS);
	//上  
	//glColor4f(1.0, 0.0, 0.0, 1.0);
	glNormal3d(0.0, 0.0, 1.0); //上  
	glVertex3d(dSize, dSize, dSize);
	glVertex3d(-dSize, dSize, dSize);
	glVertex3d(-dSize, -dSize, dSize);
	glVertex3d(dSize, -dSize, dSize);
	glEnd();
	//下  
	//glColor4f(0.0, 1.0, 0.0, 1.0);
	glBegin(GL_QUADS);
	glNormal3d(0.0, 0.0, -1.0);//下  
	glVertex3d(dSize, dSize, -dSize);
	glVertex3d(-dSize, dSize, -dSize);
	glVertex3d(-dSize, -dSize, -dSize);
	glVertex3d(dSize, -dSize, -dSize);
	glEnd();
	//前  
	//glColor4f(1.0, 1.0, 0.0, 1.0);
	glBegin(GL_QUADS);
	glNormal3d(1.0, 0.0, 0.0);//前  
	glVertex3d(dSize, dSize, dSize);
	glVertex3d(dSize, -dSize, dSize);
	glVertex3d(dSize, -dSize, -dSize);
	glVertex3d(dSize, dSize, -dSize);
	glEnd();
	//后  
	//glColor4f(0.0, 0.0, 1.0, 1.0);
	glBegin(GL_QUADS);
	glNormal3d(-1.0, 0.0, 0.0);//后  
	glVertex3d(-dSize, dSize, dSize);
	glVertex3d(-dSize, dSize, -dSize);
	glVertex3d(-dSize, -dSize, -dSize);
	glVertex3d(-dSize, -dSize, dSize);
	glEnd();
	//左  
	//glColor4f(0.0, 1.0, 0.5, 1.0);
	glBegin(GL_QUADS);
	glNormal3d(0.0, -1.0, 0.0);//左  
	glVertex3d(dSize, -dSize, dSize);
	glVertex3d(dSize, -dSize, -dSize);
	glVertex3d(-dSize, -dSize, -dSize);
	glVertex3d(-dSize, -dSize, dSize);
	glEnd();
	//右  
	//glColor4f(0.5, 1.0, 0.5, 1.0);
	glBegin(GL_QUADS);
	glNormal3d(0.0, 1.0, 0.0);//右  
	glVertex3d(dSize, dSize, dSize);
	glVertex3d(dSize, dSize, -dSize);
	glVertex3d(-dSize, dSize, -dSize);
	glVertex3d(-dSize, dSize, dSize);//右  
	glEnd();
}

/* 设置光照 */
void setLight()
{
	GLfloat light_position[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat light_ambient[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat light_diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };

	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

}

/* 绘制烟雾 */
void DrawSmoke() {

		old_smoke.CalAvgTemperature();

		for (int i = 1; i < particle_Grid_Len - 1; i++)
			for (int j = 1; j < particle_Grid_Len - 1; j++)
				for (int k = 1; k < particle_Grid_Len - 1; k++) {

					old_smoke.GetAllData(i, j, k, r, g, b,
						vx, vy, vz, vStarX, vStarY, vStarZ, size, duration, temperature, density, advectionTerm);

					double t_vox = 0;
					if (density > 0)
						t_vox = pow(e, -(extinction * h / density));
					double color = 0.8 * (1 - t_vox);
					//GLfloat mat_ambient[] = { r / 255.0, g / 255.0, b / 255.0, density};
					//GLfloat mat_diffuse[] = { r / 255.0, g / 255.0, b / 255.0, density};
					//GLfloat mat_specular[] = { r / 255.0, g / 255.0, b / 255.0, density};
					//GLfloat mat_emission[] = { 0.0f, 0.0f, 0.0f, density};
					//GLfloat mat_ambient[] = { r / 255.0, g / 255.0, b / 255.0, t_vox};
					//GLfloat mat_diffuse[] = { r / 255.0, g / 255.0, b / 255.0, t_vox};
					//GLfloat mat_specular[] = { r / 255.0, g / 255.0, b / 255.0,t_vox};
					//GLfloat mat_emission[] = { 0.0f, 0.0f, 0.0f, t_vox};

					/* 大概算是渲染？ */
					GLfloat mat_ambient[] = { r / 255.0, g / 255.0, b / 255.0, density };
					GLfloat mat_diffuse[] = { r / 255.0, g / 255.0, b / 255.0, 0.2};
					GLfloat mat_specular[] = { r / 255.0, g / 255.0, b / 255.0, density };
					GLfloat mat_emission[] = { 0.0f, 0.0f, 0.0f, density };
					GLfloat mat_shininess = 20.0f;

					//将以上材质定义应用
					glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
					glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
					glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
					glMaterialfv(GL_FRONT, GL_EMISSION, mat_emission);
					glMaterialf(GL_FRONT, GL_SHININESS, mat_shininess);

					glColor4ub(r, g, b, density * 255);
					//glColor4ub(255,255,255, density * 255);
					//glColor4f(r/255.0, g/255.0, b/255.0, 0.1f);
					//glColor4ub(r, g, b, density * 255);
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

					/* 持续时间也成影响因素*/
					if (density > 0 && duration > 0 
						&& !old_smoke.particle[i][j][k].border) {
					/* 每个烟雾粒子以球体的形式绘制 */
					glTranslated(x_location, y_location, z_location);
					glutSolidSphere(size, 10, 10);
					glTranslated(-x_location, -y_location, -z_location);
					//glFlush();
				}
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
						glEnd();*/
					//	glBegin(GL_QUADS);
					//	glVertex3d(x_location - size, y_location - size, z_location);
					//	glVertex3d(x_location - size, y_location + size, z_location);
					//	glVertex3d(x_location + size, y_location + size, z_location);
					//	glVertex3d(x_location + size, y_location - size, z_location);
					//	glEnd();
						//glFlush();
					//}
					/* 计算这个粒子所受外力 */
					old_smoke.calculateForce(i, j, k);
				}
				
		glFlush();

		for (int i = 1; i < particle_Grid_Len - 1; i++)
			for (int j = 1; j < particle_Grid_Len - 1; j++)
				for (int k = 1; k < particle_Grid_Len - 1; k++) {
					old_smoke.GetAllData(i, j, k, r, g, b,
						vx, vy, vz, vStarX,vStarY,vStarZ,size, duration, temperature, density, advectionTerm);
					/* 更新速度场 */
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
					old_smoke.SetAllData(i, j, k, r, g, b,
						vx, vy, vz, vStarX, vStarY, vStarZ, size, duration, temperature, density,advectionTerm);
					gen_smoke.SetDuration(i, j, k, duration);
				}

		semiLagrangeCalc(old_smoke, gen_smoke);

		for (int i = 1; i < particle_Grid_Len - 1; i++)
			for (int j = 1; j < particle_Grid_Len - 1; j++)
				for (int k = 1; k < particle_Grid_Len - 1; k++) {
					gen_smoke.GetAllData(i, j, k, r, g, b,
						vx, vy, vz, vStarX, vStarY, vStarZ, size, duration, temperature, density,advectionTerm);
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
					old_smoke.SetAllData(i, j, k, r, g, b,
						vx, vy, vz, vStarX, vStarY, vStarZ, size, duration - 1, temperature, density,advectionTerm);
				}

		gen_smoke.initAllData();
		glutPostRedisplay();
}

void myInit() {

	
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
	//glClearDepth(1.0f);       

	/* 开启深度测试 */
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL); 	
	//glDepthFunc(GL_ALWAYS);
	//glDepthMask(GL_FALSE);

	/* 光滑渲染 */
	glShadeModel(GL_SMOOTH);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

	/* 设置纹理 */
	glEnable(GL_TEXTURE_2D);							// Enable Texture Mapping
	glBindTexture(GL_TEXTURE_2D, texture[0]);			// Select Our Texture

	/* 混合模式 */
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	setLight();
	
	/* 两个烟雾对象的初始化 */
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

	DrawBlock(0.15);

	DrawEdge();

	DrawSmoke();

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
	//glRotatef(20.0F, 1.0F, 0.0F, 0.0F);
	//glRotatef(-30.0F, 0.0F, 1.0F, 0.0F);
	glRotatef(-20.0F, 1.0F, 0.0F, 0.0F);
	glRotatef(20.0F, 0.0F, 1.0F, 0.0F);
	//glRotatef(90.0F, 1.0F, 0.0F, 0.0F);
	//glRotatef(90.0F, 0.0F, 0.0F, 1.0F);
	//glRotatef(5.0F, 0.0F, 1.0F, 0.0F);
	//gluLookAt(0.0f, 0.5f, 0.5f, 0.5f, 0.0f, 0.0f, 0.0f, 0.0f, -0.5f);
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



