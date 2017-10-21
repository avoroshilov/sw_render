#pragma once

#include <Windows.h>

#include "types.h"
#include "math/Vec3.h"
#include "math/Mat34.h"
#include "math/Mat44.h"
#include "math/AuxMath.h"

#define DBG_USE_REGULAR_CONVERSION			0

#if (DBG_USE_REGULAR_CONVERSION == 1)
#	define m_ftoi(x) static_cast<int>(x)
#else
#	ifdef _M_AMD64
#		include <xmmintrin.h>

inline int m_ftoi(float x)
{
	return _mm_cvtt_ss2si(_mm_load_ss(&x));
}

#	else

inline int m_ftoi(float a)
{
	int retval;

	__asm fld a
	__asm fistp retval

	return retval;
}

#	endif	// _M_AMD64
#endif

// with float bilerp additional ms consumed
#define BILINEAR_FLOAT 0

class Color
{
public:

	union
	{
		struct { uchar B, G, R, A; };
		uint val;
	};

	Color()
	{
	}
	Color(const Color &c)
	{
		val = c.val;
	}

	Color operator + (const Color &c) const
	{
		Color retC;
		retC.R = R + c.R;
		retC.G = G + c.G;
		retC.B = B + c.B;
		retC.A = A + c.A;
		return retC;
	}
	Color operator - (const Color &c) const
	{
		Color retC;
		retC.R = R - c.R;
		retC.G = G - c.G;
		retC.B = B - c.B;
		retC.A = A - c.A;
		return retC;
	}

	Color & operator = (const Color &c)
	{
		val = c.val;
		return *this;
	}

	friend Color operator * (float val, const Color &c)
	{
		Color retC;
		retC.R = static_cast<uchar>(val * c.R);
		retC.G = static_cast<uchar>(val * c.G);
		retC.B = static_cast<uchar>(val * c.B);
		retC.A = static_cast<uchar>(val * c.A);
		return retC;
	}
	friend Color operator * (const Color &c, float val)
	{
		Color retC;
		retC.R = static_cast<uchar>(val * c.R);
		retC.G = static_cast<uchar>(val * c.G);
		retC.B = static_cast<uchar>(val * c.B);
		retC.A = static_cast<uchar>(val * c.A);
		return retC;
	}
};

class Texture2D
{
public:

	enum
	{
		NEAREST = 1,
		BILINEAR
	};

	uchar m_Filtration;
	uint m_Width, m_Height;
	Color *m_Bits;

	Texture2D():
		m_Bits(0),
		m_Filtration(NEAREST)
	{
	}

	~Texture2D()
	{
		if (m_Bits)
			delete [] m_Bits;
	}

	void init(uint width, uint height, uchar filtration = NEAREST)
	{
		m_Width = width;
		m_Height = height;

		m_Filtration = filtration;

		if (m_Bits)
			delete [] m_Bits;

		m_Bits = new Color[m_Width * m_Height];
	}

	Color getPoint(float x, float y)
	{
		if (x > 1.0f)
			x -= m_ftoi(x);
		if (y > 1.0f)
			y -= m_ftoi(y);

		if (x < 0.0f)
			x -= m_ftoi(x) - 1;
		if (y < 0.0f)
			y -= m_ftoi(y) - 1;

		uint ix = m_ftoi(x * (m_Width - 1)), iy = m_ftoi(y * (m_Height - 1));

		return m_Bits[ix + iy * m_Width];
	}

#if (BILINEAR_FLOAT == 1)
	Color getBilinear(float x, float y)
	{
		if (x > 1.0f)
			x -= m_ftoi(x);
		if (y > 1.0f)
			y -= m_ftoi(y);

		if (x < 0.0f)
			x -= m_ftoi(x) - 1;
		if (y < 0.0f)
			y -= m_ftoi(y) - 1;

		float fu0 = x * (m_Width - 1);
		float fv0 = y * (m_Height - 1);

		uint ix0 = m_ftoi(fu0), iy0 = m_ftoi(fv0);
		uint ix1 = (ix0 + 1) % (m_Width), iy1 = (iy0 + 1) % (m_Height);

		fu0 -= ix0;
		fv0 -= iy0;

		float fu1 = 1.0f - fu0;
		float fv1 = 1.0f - fv0;

		float W[4] =
			{
				fu1 * fv1,
				fu0 * fv1,
				fu1 * fv0,
				fu0 * fv0
			};

		Color * points[4] =
			{
				&m_Bits[ix0 + iy0 * m_Width],
				&m_Bits[ix1 + iy0 * m_Width],
				&m_Bits[ix0 + iy1 * m_Width],
				&m_Bits[ix1 + iy1 * m_Width]
			};

		Color ret;

		ret.R = points[0]->R * W[0] + points[1]->R * W[1] + points[2]->R * W[2] + points[3]->R * W[3];
		ret.G = points[0]->G * W[0] + points[1]->G * W[1] + points[2]->G * W[2] + points[3]->G * W[3];
		ret.B = points[0]->B * W[0] + points[1]->B * W[1] + points[2]->B * W[2] + points[3]->B * W[3];
		ret.A = points[0]->A * W[0] + points[1]->A * W[1] + points[2]->A * W[2] + points[3]->A * W[3];

		return ret;
	}
#else
	// Some speedboost
	Color getBilinear(float x, float y)
	{
		if (x > 1.0f)
			x -= m_ftoi(x);
		if (y > 1.0f)
			y -= m_ftoi(y);

		if (x < 0.0f)
			x -= m_ftoi(x) - 1;
		if (y < 0.0f)
			y -= m_ftoi(y) - 1;

		float mappedCoordX = x * (m_Width - 1);
		float mappedCoordY = y * (m_Height - 1);

		uint ix0 = m_ftoi(mappedCoordX), iy0 = m_ftoi(mappedCoordY);
		uint ix1 = (ix0 + 1) % (m_Width), iy1 = (iy0 + 1) % (m_Height);

		uchar fu0 = static_cast<uchar>((mappedCoordX - ix0) * 255);
		uchar fv0 = static_cast<uchar>((mappedCoordY - iy0) * 255);

		uint fu1 = 255 - fu0;
		uint fv1 = 255 - fv0;

		iy0 *= m_Width;
		iy1 *= m_Width;

		uint W[4] =
			{
				((uint)fu1 * fv1) >> 8,
				((uint)fu0 * fv1) >> 8,
				((uint)fu1 * fv0) >> 8,
				((uint)fu0 * fv0) >> 8
			};

		Color * points[4] =
			{
				&m_Bits[ix0 + iy0],
				&m_Bits[ix1 + iy0],
				&m_Bits[ix0 + iy1],
				&m_Bits[ix1 + iy1]
			};

		Color ret;

		ret.R = (points[0]->R * W[0] + points[1]->R * W[1] + points[2]->R * W[2] + points[3]->R * W[3]) >> 8;
		ret.G = (points[0]->G * W[0] + points[1]->G * W[1] + points[2]->G * W[2] + points[3]->G * W[3]) >> 8;
		ret.B = (points[0]->B * W[0] + points[1]->B * W[1] + points[2]->B * W[2] + points[3]->B * W[3]) >> 8;
		ret.A = (points[0]->A * W[0] + points[1]->A * W[1] + points[2]->A * W[2] + points[3]->A * W[3]) >> 8;

		return ret;
	}
#endif

	Color get(float x, float y)
	{
		if (m_Filtration == NEAREST)
		{
			return getPoint(x, y);
		}
		else
		{
			return getBilinear(x, y);
		}
	}

};

class Vertex
{
public:

	float inv_w;
	math::Vec4 fColor;

	float w;

	math::Vec3 coord;
	math::Vec3 normal;
	math::Vec3 texCoord;
	Color color;

	Vertex()
	{
	}

	Vertex(const Vertex &v)
	{
		coord = v.coord;
		normal = v.normal;
		texCoord = v.texCoord;
		color = v.color;
		fColor = v.fColor;
		w = v.w;
		inv_w = v.inv_w;
	}

	Vertex & operator = (const Vertex &v)
	{
		coord = v.coord;
		normal = v.normal;
		texCoord = v.texCoord;
		color = v.color;
		fColor = v.fColor;
		w = v.w;
		inv_w = v.inv_w;

		return *this;
	}
};

class Drawer
{
	uint tsCount;
	math::Mat34 * transformStack;
	uint psCount;
	math::Mat44 * projectionStack;

public:

	math::Mat34 transform;
	math::Mat44 projection;

	Color color;
	math::Vec3 normal;
	math::Vec3 texCoord;

	bool useSphericalEnvMapping;

	Texture2D *tex2D;

	uint m_width, m_height;

	uint *m_bits;
	uint *m_zBuffer;

	HDC m_hDC, m_hExtDC;
	HBITMAP m_hBitmap;
	HBITMAP m_hInitialBitmap;

	Drawer():
		m_zBuffer(0),
		m_bits(0),
		transformStack(0),
		projectionStack(0),
		tsCount(0),
		psCount(0),
		tex2D(0),
		useSphericalEnvMapping(false)
	{
		transformStack = new math::Mat34[64];
		projectionStack = new math::Mat44[64];
	}

	~Drawer()
	{
		deinit();

		if (transformStack)
			delete [] transformStack;
		if (projectionStack)
			delete [] projectionStack;
	}

	void pushTransform()
	{
		transformStack[tsCount] = transform;
		++tsCount;
	}
	void pullTransform()
	{
		--tsCount;
		transform = transformStack[tsCount];
	}

	void pushProjection()
	{
		projectionStack[psCount] = projection;
		++psCount;
	}
	void pullProjection()
	{
		--psCount;
		projection = projectionStack[psCount];
	}

	math::Mat34 & getTransformMat()
	{
		return transform;
	}
	math::Mat44 & getProjectionMat()
	{
		return projection;
	}

	void projPerspective(float FOV, float asp, float N, float F, float W = 1.0f, float H = 1.0f)
	{
		projection.zero();

		W = 1.0f / (N * tanf(FOV));
		H = asp * W;

		projection._00 = N * W;
		projection._11 = N * H;
		projection._22 = -1.0f / (F - N);
		projection._23 = -N / (F - N);
		projection._32 = -1.0f;
	}

	void resetBuffers()
	{
		memset(m_bits, 0, m_width * m_height * 4 * sizeof(uchar));
		memset(m_zBuffer, 0xFF, m_width * m_height * sizeof(uint));
	}

	void swapBuffers()
	{
		BitBlt(
			m_hExtDC,
			0,
			0,
			m_width, 
			m_height,
			m_hDC,
			0,
			0,
			SRCCOPY
			);
	}

	void renderDepth()
	{
		for (uint i = 0; i < m_width; ++i)
		{
			for (uint j = 0; j < m_height; ++j)
			{
				uint val = uint((m_zBuffer[i + j*m_width]-1) / 16843009.0f);
				m_bits[i + j*m_width] = (val << 16) + (val << 8) + val;
			}
		}
	}

	void setColorF(float R, float G, float B, float A = 1.0f)
	{
		color.R = static_cast<uchar>(255.0f * R);
		color.G = static_cast<uchar>(255.0f * G);
		color.B = static_cast<uchar>(255.0f * B);
		color.A = static_cast<uchar>(255.0f * A);
	}
	void setColorUC(uchar R, uchar G, uchar B, uchar A = 255)
	{
		color.R = R;
		color.G = G;
		color.B = B;
		color.A = A;
	}

	math::Vec3 toViewportCopy(math::Vec3 modelCoords)
	{
		// to world
		transform.transform(modelCoords);
		// to clip space
		projection.homogTransform(modelCoords);

		// to window coords [viewport]
		modelCoords.x += 1.0f;
		modelCoords.y += 1.0f;
		modelCoords.x *= 0.5f * m_width;
		modelCoords.y *= 0.5f * m_height;

		return modelCoords;
	}

	math::Vec3 & toViewport(math::Vec3 &modelCoords)
	{
		// to world
		transform.transform(modelCoords);
		// to clip space
		projection.homogTransform(modelCoords);

		// to window coords [viewport]
		modelCoords.x += 1.0f;
		modelCoords.y += 1.0f;
		modelCoords.x *= 0.5f * m_width;
		modelCoords.y *= 0.5f * m_height;

		return modelCoords;
	}

	Vertex toVertex(math::Vec3 modelCoords, uint flags = 0);
	Vertex toVertexLoc(math::Vec3 modelCoords, uint flags = 0);

	void putPixel_us(int x, int y, uint color)
	{
		m_bits[x + y * m_width] = color;
	}

	void drawLineFogged(const Vertex &p0, const Vertex &p1);

	// Function interpolates vertex attributes and performs scanline plot (X axis)
	// Suffixes:
	//	*R - attribute root (all the interpolation is done relative to this element)
	//	*R1 - diff in attribute between root point and point 1
	//	*RO - diff in attribute between other remaining root point and current root point
	//			i.e. in case root point is 0 - remaing is 2; in case root is 2 - remainnig is 0
	// Attributes:
	//	sp* - sorted point, basically just vertex coordinate
	//	inv_w* - 1.0/w (perspective division) - only for rendering with perspective correction
	//	c* - color
	//	n* - normal
	//	tc* - texture coordinates
	void interpolateAndPlotScanline(
		int y,

		const math::Vec3 &spR,
		const math::Vec3 &dR1,
		const math::Vec3 &dRO,

		const float *inv_wR,
		float inv_wR1,
		float inv_wRO,

		const math::Vec4 *cR,
		const math::Vec4 &cR1,
		const math::Vec4 &cRO,

		const math::Vec3 *nR,
		const math::Vec3 &nR1,
		const math::Vec3 &nRO,

		const math::Vec3 *tcR,
		const math::Vec3 &tcR1,
		const math::Vec3 &tcRO
		);
	void drawTriangle(const Vertex &p0, const Vertex &p1, const Vertex &p2);
	
	void renderTriangle(Vertex p0, Vertex p1, Vertex p2);

	void setParams(uint width, uint height)
	{
		m_width = width;
		m_height = height;
	}

	void init(HDC extHDC)
	{
		m_hExtDC = extHDC;

		BITMAPINFO *bmInfo = new BITMAPINFO;
		memset(bmInfo, 0, sizeof(BITMAPINFO));

		bmInfo->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
		bmInfo->bmiHeader.biWidth = m_width;
		bmInfo->bmiHeader.biHeight = m_height;
		bmInfo->bmiHeader.biPlanes = 1;
		bmInfo->bmiHeader.biBitCount = 32;
		bmInfo->bmiHeader.biCompression = BI_RGB;
		bmInfo->bmiHeader.biSizeImage = m_width * m_height * 4;

		if (m_zBuffer != 0)
			delete [] m_zBuffer;
		
		m_zBuffer = new uint[m_width * m_height];

		m_hBitmap = CreateDIBSection(m_hExtDC, bmInfo, DIB_RGB_COLORS, (void **)&m_bits, NULL, 0x0);
		m_hDC = CreateCompatibleDC(m_hExtDC);
		m_hInitialBitmap = (HBITMAP)SelectObject(m_hDC, m_hBitmap);
	}

	void deinit()
	{
		if (m_hBitmap)
		{
			if (m_hDC)
			{
				if (m_hInitialBitmap)
				{
					SelectObject(m_hDC, m_hInitialBitmap);
					m_hInitialBitmap = 0;
				}

				DeleteDC(m_hDC);
				m_hDC = 0;
			}

			DeleteObject(m_hBitmap);
			m_hBitmap = 0;
		}

		m_bits = 0;
		if (m_zBuffer != 0)
		{
			delete [] m_zBuffer;
			m_zBuffer = 0;
		}

		m_width = 0;
		m_height = 0;
	}
};
