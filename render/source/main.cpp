#include <windows.h>

#include <stdio.h>

#include "types.h"
#include "windows/window.h"
#include "windows/timer.h"
#include "math/Vec3.h"
#include "math/mat34.h"
#include "math/mat44.h"
#include "drawer.h"

using namespace math;

bool isRunning = true;

windows::Window mainWindow;
Drawer mainDrawer;

Texture2D checker;

float Time = 0.0f;

bool isPaused = false;

const uint numSegm = 64;
const uint numVertices = numSegm * numSegm;
Vec3 spVerts[numVertices], spNormals[numVertices];

void Render()
{
	mainDrawer.resetBuffers();

	float radius = 0.7f;
	float amplitude = 0.2f, phase = Time / 500.0f;
	Vec3 Sines = Vec3C(2, 2, 3);
	Vec3 offset = Vec3C(0.0f, 0.0f, 0.0f);

	Vec3 tmp;
	for (uint j = 0; j < numSegm; j++)
	{
		float AngJ = j*PI/(float)(numSegm - 1) - _PI2;

		for (int i = 0; i < numSegm; i++)
		{
			float AngI = i*_2PI/(float)(numSegm-1);

			Vec3 *spPoint = &spVerts[i + j*numSegm];

			*spPoint = offset + Vec3C(radius * cosf(AngJ)*cosf(AngI), radius * sinf(AngJ), radius * cosf(AngJ)*sinf(AngI));

			tmp.x = sinf(spPoint->y*spPoint->z * PI * Sines.x + phase) * amplitude;
			tmp.y = cosf(spPoint->x*spPoint->z * PI * Sines.y + phase) * amplitude;
			tmp.z = sinf(spPoint->x*spPoint->y * PI * Sines.z + phase) * amplitude;

			*spPoint += tmp;
		}
	}
	for (uint cnt = 0; cnt < numSegm*numSegm; cnt++)
	{
		spNormals[cnt] = Vec3C(0.0f, 0.0f, 0.0f);
	}
	for (uint j = 0; j < numSegm; j++)
	{
		uint jAdd0 = j * numSegm;
		uint jAdd1 = ((j + 1) % numSegm) * numSegm;

		for (int i0 = 0; i0 < numSegm; ++i0)
		{
			uint i1 = i0 + 1;
			if (i1 == numSegm)
			{
				// numSegm-1 and 0 vertices exactly match, so diff would be 0 and normal incorrect, need to use slice#1
				i1 = 1;
			}
			Vec3 normal = (spVerts[i0+jAdd1] - spVerts[i0+jAdd0]).cross(spVerts[i1+jAdd0] - spVerts[i0+jAdd0]).getNormalized();
			spNormals[i0 + jAdd0] += normal;
		}
	}
	for (uint cnt = 0; cnt < numSegm*numSegm; cnt++)
	{
		spNormals[cnt].normalize();
	}

	Mat34 &transform = mainDrawer.getTransformMat();
	transform.identity();
	transform.fillTranslation(Vec3C(0.0f, 0.0f, -3.0f));

	mainDrawer.pushTransform();

	transform.rotate(Vec3C(0, 1, 0), Time / 1000.0f);

	Vertex v0, v1, v2, v_tmp0, v_tmp1;

	mainDrawer.tex2D = &checker;
	mainDrawer.useSphericalEnvMapping = false;

#define RENDER_DEFORM_SPHERE	1
#define RENDER_WIRED_CUBE		1
#define RENDER_SOLID_CUBE		0

#if (RENDER_DEFORM_SPHERE == 1)

#define RENDER_DEFORM_SPHERE_NORMALS	1

	mainDrawer.setColorUC(255, 255, 255);
	mainDrawer.texCoord = Vec3C(0.0f, 0.0f, 0.0f);

	// 2k tris
	for (uint j = 0; j < numSegm - 1; j++)
	{
		uint jAdd0 = j * numSegm;
		uint jAdd1 = (j + 1) * numSegm;

		float tcJ0 = (j  ) / (float)numSegm;
		float tcJ1 = (j+1) / (float)numSegm;

		for (int i0 = 0; i0 < numSegm - 1; ++i0)
		{
			uint i1 = (i0 + 1) % numSegm;

			//mainDrawer.setColorUC(255, 0, 0);
			mainDrawer.normal = spNormals[i0+jAdd0];
			mainDrawer.texCoord.x = i0 / (float)numSegm;
			mainDrawer.texCoord.y = tcJ0;
			v0 = mainDrawer.toVertexLoc(spVerts[i0+jAdd0]);

			//mainDrawer.setColorUC(0, 255, 0);
			mainDrawer.normal = spNormals[i1+jAdd0];
			mainDrawer.texCoord.x = (i0+1) / (float)numSegm;
			mainDrawer.texCoord.y = tcJ0;
			v1 = mainDrawer.toVertexLoc(spVerts[i1+jAdd0]);

			//mainDrawer.setColorUC(0, 0, 255);
			mainDrawer.normal = spNormals[i0+jAdd1];
			mainDrawer.texCoord.x = i0 / (float)numSegm;
			mainDrawer.texCoord.y = tcJ1;
			v2 = mainDrawer.toVertexLoc(spVerts[i0+jAdd1]);

			mainDrawer.renderTriangle(v0, v1, v2);

#if (RENDER_DEFORM_SPHERE_NORMALS == 1)
			v_tmp0 = mainDrawer.toVertex(spVerts[i0+jAdd0]);
			v_tmp1 = mainDrawer.toVertex(spVerts[i0+jAdd0] + 0.1f * mainDrawer.normal);
			mainDrawer.drawLineFogged(v_tmp0, v_tmp1);
#endif

			mainDrawer.normal = spNormals[i1+jAdd1];

			//mainDrawer.setColorUC(255, 0, 0);
			mainDrawer.texCoord.x = (i0+1) / (float)numSegm;
			mainDrawer.texCoord.y = tcJ1;
			v0 = mainDrawer.toVertexLoc(spVerts[i1+jAdd1]);

			mainDrawer.renderTriangle(v0, v2, v1);
		}
	}

#endif

	mainDrawer.pullTransform();

#if (RENDER_WIRED_CUBE == 1)

	mainDrawer.pushTransform();

	transform.fillRotation(Vec3C(0.7f, 0.5f, 1.0f).getNormalized(), phase * 0.1f);
	transform.rotate(Vec3C(1.0f, 0.8f, 0.2f).getNormalized(), phase * 0.3f);

	float cubeSize = 0.5f;

	Vec3 cubeCoords[8] =
		{
			Vec3C(-cubeSize, -cubeSize, -cubeSize),	// LLL [0]
			Vec3C( cubeSize, -cubeSize, -cubeSize),	// ULL [1]
			Vec3C(-cubeSize,  cubeSize, -cubeSize),	// LUL [2]
			Vec3C( cubeSize,  cubeSize, -cubeSize),	// UUL [3]
			Vec3C(-cubeSize, -cubeSize,  cubeSize),	// LLU [4]
			Vec3C( cubeSize, -cubeSize,  cubeSize),	// ULU [5]
			Vec3C(-cubeSize,  cubeSize,  cubeSize),	// LUU [6]
			Vec3C( cubeSize,  cubeSize,  cubeSize),	// UUU [7]
		};

	uchar cubeColors[8][3] = 
		{
			{ 255, 0, 0},
			{ 0, 255, 0},
			{ 0, 0, 255},
			{ 255, 128, 0},
			{ 255, 255, 0},
			{ 0, 128, 255},
			{ 0, 255, 255},
			{ 255, 255, 255},
		};

	Vertex cubeVerts[8];

	for (uint i = 0; i < 8; ++i)
	{
		mainDrawer.setColorUC(cubeColors[i][0], cubeColors[i][1], cubeColors[i][2]);
		cubeVerts[i] = mainDrawer.toVertex(cubeCoords[i]);
	}

#define CUBE_DRAW_LINE(idx0, idx1)\
		mainDrawer.drawLineFogged(cubeVerts[idx0], cubeVerts[idx1]);

	// Back face
	CUBE_DRAW_LINE(0, 1);
	CUBE_DRAW_LINE(1, 3);
	CUBE_DRAW_LINE(3, 2);
	CUBE_DRAW_LINE(2, 0);
	// Front face
	CUBE_DRAW_LINE(4, 5);
	CUBE_DRAW_LINE(5, 7);
	CUBE_DRAW_LINE(7, 6);
	CUBE_DRAW_LINE(6, 4);
	// Left
	CUBE_DRAW_LINE(6, 2);
	CUBE_DRAW_LINE(4, 0);
	// Right
	CUBE_DRAW_LINE(7, 3);
	CUBE_DRAW_LINE(5, 1);

#undef CUBE_DRAW_LINE

	mainDrawer.pullTransform();

#endif


#if (RENDER_SOLID_CUBE == 1)
	mainDrawer.pushTransform();

	mainDrawer.setColorUC(255, 255, 255, 255);

	transform.fillRotation(Vec3C(0.7f, 0.5f, 1.0f).getNormalized(), phase * 0.1f);
	transform.rotate(Vec3C(1.0f, 0.8f, 0.2f).getNormalized(), phase * 0.3f);

	Drawer & refDrawer = mainDrawer;
	auto drawFace = [&refDrawer](const Vec3 & normal, float size1, float size2, float sizeN)
	{
		Vertex v0, v1, v2;

		Vec3 ax1, ax2;
		normal.tangentSpace(ax1, ax2);

#if 0
#	define SET_COLOR(r, g, b) 	mainDrawer.setColorUC(r, g, b);
#else
#	define SET_COLOR(r, g, b)
#endif

		mainDrawer.normal = normal;
		SET_COLOR(255, 255, 255);
		mainDrawer.texCoord.x = 0.0f;
		mainDrawer.texCoord.y = 0.0f;
		v0 = mainDrawer.toVertexLoc(normal*sizeN - ax1*size1 - ax2*size2);
		SET_COLOR(255, 0, 255);
		mainDrawer.texCoord.x = 0.0f;
		mainDrawer.texCoord.y = 1.0f;
		v1 = mainDrawer.toVertexLoc(normal*sizeN - ax1*size1 + ax2*size2);
		SET_COLOR(0, 0, 255);
		mainDrawer.texCoord.x = 1.0f;
		mainDrawer.texCoord.y = 1.0f;
		v2 = mainDrawer.toVertexLoc(normal*sizeN + ax1*size1 + ax2*size2);
		mainDrawer.renderTriangle(v0, v1, v2);

		SET_COLOR(255, 255, 255);
		mainDrawer.texCoord.x = 0.0f;
		mainDrawer.texCoord.y = 0.0f;
		v0 = mainDrawer.toVertexLoc(normal*sizeN - ax1*size1 - ax2*size2);
		SET_COLOR(0, 0, 255);
		mainDrawer.texCoord.x = 1.0f;
		mainDrawer.texCoord.y = 1.0f;
		v1 = mainDrawer.toVertexLoc(normal*sizeN + ax1*size1 + ax2*size2);
		SET_COLOR(0, 255, 255);
		mainDrawer.texCoord.x = 1.0f;
		mainDrawer.texCoord.y = 0.0f;
		v2 = mainDrawer.toVertexLoc(normal*sizeN + ax1*size1 - ax2*size2);
		mainDrawer.renderTriangle(v0, v1, v2);
	};
	
	const Vec3 cubeSizeVec = Vec3C(0.7f, 0.7f, 0.7f);
	drawFace(Vec3C( 0.0f, 0.0f, 1.0f), cubeSizeVec.y, cubeSizeVec.x, cubeSizeVec.z);
	drawFace(Vec3C( 0.0f, 0.0f,-1.0f), cubeSizeVec.y, cubeSizeVec.x, cubeSizeVec.z);
	drawFace(Vec3C( 0.0f, 1.0f, 0.0f), cubeSizeVec.x, cubeSizeVec.z, cubeSizeVec.y);
	drawFace(Vec3C( 0.0f,-1.0f, 0.0f), cubeSizeVec.x, cubeSizeVec.z, cubeSizeVec.y);
	drawFace(Vec3C( 1.0f, 0.0f, 0.0f), cubeSizeVec.y, cubeSizeVec.z, cubeSizeVec.x);
	drawFace(Vec3C(-1.0f, 0.0f, 0.0f), cubeSizeVec.y, cubeSizeVec.z, cubeSizeVec.x);

	mainDrawer.pullTransform();
#endif

	mainDrawer.tex2D = 0;

	//mainDrawer.renderDepth();
	mainDrawer.swapBuffers();
}

void BasicInit()
{
	uint width, height;
	width = mainWindow.getWidth();
	height = mainWindow.getHeight();
	mainDrawer.projPerspective(PI / 4.0f, width / (float)height, 0.01f, 10.0f, 0.1f, 0.1f);

	const int xTS = 256;
	const int yTS = 256;
	unsigned char * image = new uchar [xTS*yTS*4];
	unsigned char * temp = new uchar [xTS*yTS*4];

	uchar some_val = 2;
	{
		int yAdd = 0;
		for (int j = 0; j < yTS; ++j)
		{
			for (int i = 0; i < xTS; ++i)
			{
				uchar c = ((i << some_val) ^ (j << some_val)) >> 2;
				image[(i+yAdd)*4+0] = c;
				image[(i+yAdd)*4+1] = c;
				image[(i+yAdd)*4+2] = c;
				image[(i+yAdd)*4+3] = c;
			}
			yAdd += xTS;
		}
	}

	if (1)
	{
		int x, y;
		float R, an;
		int mTS = xTS >> 1;

		short lX = (xTS >> 1), lY = (yTS >> 1);

		float lE = 1.1f;
		float sqrY, subY;

		float PI_2 = PI / 2.0f;

		unsigned char *lPos = image, *tPos = temp;
		for (y = 0; y < yTS; y++)
		{
			sqrY = float(sqr(y - lY));
			subY = float((yTS-y) - (lY));
			for (x = 0; x < xTS; x++)
			{
				R = sqrtf(float(sqr(x - lX)) + sqrY);
				an = atan2f(subY, float(x - (lX))) - PI_2;
				if (R < mTS) R = mTS * expf(logf(R / mTS) * lE);

				//GetBerpTexel_RGB(lPos, lX + (R*cosf(an)), lY + (R*sinf(an)), tPos, tPos + 1, tPos + 2);
				int xCoord = (int)(lX + (R*cosf(an)));
				int yCoord = (int)(lY + (R*sinf(an)));
				xCoord = xCoord % xTS;
				yCoord = yCoord % yTS;

				*(tPos+0) = *(lPos+(xCoord+yCoord*xTS)*4+0);
				*(tPos+1) = *(lPos+(xCoord+yCoord*xTS)*4+1);
				*(tPos+2) = *(lPos+(xCoord+yCoord*xTS)*4+2);
				tPos += 4;
			}
		}

		memcpy(image, temp, (xTS*yTS<<2)*sizeof(unsigned char));
	}

	checker.init(xTS, yTS, Texture2D::BILINEAR);
	{
		int yAdd = 0;
		for (int j = 0; j < yTS; ++j)
		{
			for (int i = 0; i < xTS; ++i)
			{
				Color color;
				color.R = image[(i+yAdd)*4+0];
				color.G = image[(i+yAdd)*4+1];
				color.B = image[(i+yAdd)*4+2];
				color.A = image[(i+yAdd)*4+3];
				checker.m_Bits[i + j * xTS] = color;
			}
			yAdd += xTS;
		}
	}

	delete [] image;
	delete [] temp;
}

void BasicDeinit()
{
}

void chageFocusCallback(void * pUserData, bool isInFocus)
{
}
void keyStateCallback(void * pUserData, windows::KeyCode keyCode, windows::KeyState keyState)
{
	using namespace windows;

	switch (keyCode)
	{
		case KeyCode::eEscape:
		{
			isRunning = false;
			break;
		}
		case (KeyCode)'P':
		{
			if (keyState == KeyState::ePressed)
			{
				isPaused = !isPaused;
			}
			//processed = true;
			break;
		}
		case (KeyCode)'O':
		{
			if (keyState == KeyState::ePressed || keyState == KeyState::eHeldDown)
			{
				isPaused = true;
				Time += 16.0f;
			}
			//processed = true;
		}
	}
}

int WINAPI WinMain(HINSTANCE hInst,HINSTANCE hPrevInst,
				   LPSTR cmdLine,int cmdShow)
{
	using namespace windows;

	uint width = 640, height = 480, bpp = 32;

	mainWindow.setResizeable(false);
	mainWindow.setParameters(width, height, Window::Kind::eWindowed);
	mainWindow.setTitle(L"Window title");
	mainWindow.init(true);

	mainDrawer.setParams(width, height);
	mainDrawer.init(mainWindow.getDC());

	setChangeFocusCallback(chageFocusCallback);
	setKeyStateCallback(keyStateCallback);

	BasicInit();

	isRunning = true;

	windows::Timer fpsTimer;
	fpsTimer.start();

	uint seconds = 0;

	float dts = 0.0f;
	uint frames = 0;

	MSG msg;
	while (isRunning)
	{
		if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
		{
			if (msg.message == WM_QUIT)
				break;

			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}

//		mainDrawer.test();
		Render();

		double dt = fpsTimer.time();
		fpsTimer.start();

		if (!isPaused)
		{
			Time += (float)dt;
		}

		dts += (float)dt;
		++frames;

		if (dts > 500.0f)
		{
			float av_dt = dts / frames;
			
			dts = 0.0f;
			frames = 0;

			float FPS = float(1000.0f / av_dt);
			mainWindow.setTitle(L"Window title [%d] - %.2fms", (int)FPS, dt);
		}
	}

	BasicDeinit();
	mainWindow.deinit();

	return 0;
}
