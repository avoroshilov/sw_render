#include "drawer.h"

#define PERSPECTIVE_CORRECTION 1

// 0 - no normalization		[0]
// 1 - fast normalization	[1e-2f]
// 2 - real normalization	[1e-5f]
#define NORMALIZATION 1
//#define NORMALIZATION_EPS 1e-2f

#define BACKFACE_CULL_CW				0
#define BACKFACE_CULL_CCW				1	// default

// 0 - usual color
// 1 - normal
// 2 - texture coordinate
#define DBG_TRIANGLE_RASTER				0
#define DBG_RASTERIZE_EDGES				0

using namespace math;

void SutherlandHodgman(const Vertex & v0, const Vertex & v1, bool ok0, bool ok1, Vertex * clipped, uint & clippedVertices, bool isLastEdge = false)
{
	if (ok0 && ok1)
	{
		clipped[clippedVertices++] = v0;
		
		//if (!isLastEdge)
		//	clipped[clippedVertices++] = v1;
	}
	else if (ok0 && !ok1)
	{
		clipped[clippedVertices++] = v0;
		
		// Intersection 0-1 with plane z = 0
		Vec3 diff = v1.coord - v0.coord;
		float t = -v0.coord.z / diff.z;

		Vertex & refVertex = clipped[clippedVertices++];

		refVertex.coord = v0.coord + diff * t;
		refVertex.normal = v0.normal + (v1.normal - v0.normal) * t;
		refVertex.texCoord = v0.texCoord + (v1.texCoord - v0.texCoord) * t;
		refVertex.w = v0.w + (v1.w - v0.w) * t;
	
		refVertex.color.val = v0.color.val;
		refVertex.color.R += m_ftoi(v1.color.R - v0.color.R * t);
		refVertex.color.G += m_ftoi(v1.color.G - v0.color.G * t);
		refVertex.color.B += m_ftoi(v1.color.B - v0.color.B * t);
		refVertex.color.A += m_ftoi(v1.color.A - v0.color.A * t);
	}
	else if (!ok0 && ok1)
	{
		// Intersection 0-1 with plane z = 0
		Vec3 diff = v1.coord - v0.coord;
		float t = -v0.coord.z / diff.z;

		Vertex & refVertex = clipped[clippedVertices++];

		refVertex.coord = v0.coord + diff * t;
		refVertex.normal = v0.normal + (v1.normal - v0.normal) * t;
		refVertex.texCoord = v0.texCoord + (v1.texCoord - v0.texCoord) * t;
		refVertex.w = v0.w + (v1.w - v0.w) * t;
	
		refVertex.color.val = v0.color.val;
		refVertex.color.R += m_ftoi(v1.color.R - v0.color.R * t);
		refVertex.color.G += m_ftoi(v1.color.G - v0.color.G * t);
		refVertex.color.B += m_ftoi(v1.color.B - v0.color.B * t);
		refVertex.color.A += m_ftoi(v1.color.A - v0.color.A * t);
	}
}

Vertex Drawer::toVertex(math::Vec3 modelCoords, uint flags)
{
	Vertex ret;
	ret.coord = modelCoords;

	// to world
	transform.transform(ret.coord);
	// to clip space
	projection.homogTransform(ret.coord, ret.w);

	// perspective division
	float inv_w = 1.0f / ret.w;
	ret.inv_w = inv_w;
	ret.coord.x *= inv_w;
	ret.coord.y *= inv_w;
	ret.coord.z *= inv_w;

	// to window coords [viewport]
	ret.coord.x += 1.0f;
	ret.coord.y += 1.0f;
	ret.coord.x *= 0.5f * m_width;
	ret.coord.y *= 0.5f * m_height;

	math::Vec3 trNormal(normal);
	// Transform normal
	//	n^T*v = 0
	//	n^T*(M^-1)*M*v = 0
	//	=> n'^T = n^T*(M^-1)
	//	n' = (M^-T)*n
	// Since rotation matrix is orthogonal
	ret.normal = transform.rotate(trNormal);
	ret.normal.normalize();

	memcpy((void *)&ret.color, (void *)&color, sizeof(Color));
	ret.texCoord = texCoord;

#if (PERSPECTIVE_CORRECTION == 1)
	// prepare perspective correction
	ret.fColor.x = (color.R*inv_w);
	ret.fColor.y = (color.G*inv_w);
	ret.fColor.z = (color.B*inv_w);
	ret.fColor.w = (color.A*inv_w);
	ret.normal *= inv_w;
	ret.texCoord *= inv_w;
#else
	ret.fColor.x = (float)color.R;
	ret.fColor.y = (float)color.G;
	ret.fColor.z = (float)color.B;
	ret.fColor.w = (float)color.A;
#endif

	// Converting to 0..1 for convenience
	ret.fColor.x /= 255.0f;
	ret.fColor.y /= 255.0f;
	ret.fColor.z /= 255.0f;
	ret.fColor.w /= 255.0f;

	return ret;
}

Vertex Drawer::toVertexLoc(math::Vec3 modelCoords, uint flags)
{
	Vertex ret;
	ret.coord = modelCoords;
	ret.normal = normal;

	memcpy((void *)&ret.color, (void *)&color, sizeof(Color));
	ret.texCoord = texCoord;

	return ret;
}

void Drawer::renderTriangle(Vertex v0, Vertex v1, Vertex v2)
{
	Mat44 MVP = projection * transform;

	// to world, to clip space
	MVP.homogTransform(v0.coord, v0.w);
	MVP.homogTransform(v1.coord, v1.w);
	MVP.homogTransform(v2.coord, v2.w);

	// Transform normal
	//	n^T*v = 0
	//	n^T*(M^-1)*M*v = 0
	//	=> n'^T = n^T*(M^-1)
	//	n' = (M^-T)*n
	// Since rotation matrix is orthogonal
	transform.rotate(v0.normal);
	transform.rotate(v1.normal);
	transform.rotate(v2.normal);

	// Sutherland-Hodgman clipping
	uint clippedVertices = 0;
	Vertex clipped[4];
	
	bool isOk[3];
	isOk[0] = (v0.coord.z > 0.0f);
	isOk[1] = (v1.coord.z > 0.0f);
	isOk[2] = (v2.coord.z > 0.0f);

	// 0-1
	SutherlandHodgman(v0, v1, isOk[0], isOk[1], clipped, clippedVertices, false);
	// 1-2
	SutherlandHodgman(v1, v2, isOk[1], isOk[2], clipped, clippedVertices, false);
	// 2-0
	SutherlandHodgman(v2, v0, isOk[2], isOk[0], clipped, clippedVertices, true);
	
	float hW = 0.5f * m_width, hH = 0.5f * m_height;

	// perspective division
	for (uint i = 0; i < clippedVertices; ++i)
	{
		Vertex & refVertex = clipped[i];

		refVertex.inv_w = 1.0f / refVertex.w;
		refVertex.coord.x *= refVertex.inv_w;
		refVertex.coord.y *= refVertex.inv_w;
		refVertex.coord.z *= refVertex.inv_w;

		// to viewport
		refVertex.coord.x += 1.0f;
		refVertex.coord.y += 1.0f;
		refVertex.coord.x *= hW;
		refVertex.coord.y *= hH;

#if (PERSPECTIVE_CORRECTION == 1)
		// prepare perspective correction
		refVertex.fColor.x = (refVertex.color.R*refVertex.inv_w);
		refVertex.fColor.y = (refVertex.color.G*refVertex.inv_w);
		refVertex.fColor.z = (refVertex.color.B*refVertex.inv_w);
		refVertex.fColor.w = (refVertex.color.A*refVertex.inv_w);
		refVertex.normal *= refVertex.inv_w;
		refVertex.texCoord *= refVertex.inv_w;
#else
		refVertex.fColor.x = (float)refVertex.color.R;
		refVertex.fColor.y = (float)refVertex.color.G;
		refVertex.fColor.z = (float)refVertex.color.B;
		refVertex.fColor.w = (float)refVertex.color.A;
#endif

		// Converting to 0..1 for convenience
		refVertex.fColor.x /= 255.0f;
		refVertex.fColor.y /= 255.0f;
		refVertex.fColor.z /= 255.0f;
		refVertex.fColor.w /= 255.0f;
	}

	if (clippedVertices == 3)
	{
		// Only one triangle needs to be rasterized
		drawTriangle(clipped[0], clipped[1], clipped[2]);
	}
	else if (clippedVertices == 4)
	{
		// Cut triangle forms a quad, need to be split into two triangles
		drawTriangle(clipped[0], clipped[1], clipped[3]);
		drawTriangle(clipped[1], clipped[2], clipped[3]);
	}
}

void Drawer::drawLineFogged(const Vertex &v0, const Vertex &v1)
{
	int temp;

	int x0 = m_ftoi(v0.coord.x), y0 = m_ftoi(v0.coord.y);
	int x1 = m_ftoi(v1.coord.x), y1 = m_ftoi(v1.coord.y);
	float z0 = v0.coord.z, z1 = v1.coord.z;

	int deltax = x1 - x0;
	int deltay = y1 - y0;

	int pDeltaX = deltax;
	if (deltax < 0)
		pDeltaX = -deltax;
	int pDeltaY = deltay;
	if (deltay < 0)
		pDeltaY = -deltay;

	float deltaz = z1 - z0;

#if (PERSPECTIVE_CORRECTION == 1)
	float inv_w0 = v0.inv_w;
	float delta_inv_w = v1.inv_w - v0.inv_w;
#endif

	int deltaC[4];
	const Color *c0;

	if (pDeltaX > pDeltaY)
	{
		if (x0 > x1)
		{
			temp = x0;
			x0 = x1;
			x1 = temp;

			temp = y0;
			y0 = y1;
			y1 = temp;

			float tempF = z0;
			z0 = z1;
			z1 = tempF;

			deltax = -deltax;
			deltay = -deltay;
			deltaz = -deltaz;

#if (PERSPECTIVE_CORRECTION == 1)
			inv_w0 = v1.inv_w;
			delta_inv_w = -delta_inv_w;
#endif

			deltaC[0] = v0.color.R - v1.color.R;
			deltaC[1] = v0.color.G - v1.color.G;
			deltaC[2] = v0.color.B - v1.color.B;
			deltaC[3] = v0.color.A - v1.color.A;
			c0 = &v1.color;
		}
		else
		{
			deltaC[0] = v1.color.R - v0.color.R;
			deltaC[1] = v1.color.G - v0.color.G;
			deltaC[2] = v1.color.B - v0.color.B;
			deltaC[3] = v1.color.A - v0.color.A;
			c0 = &v0.color;
		}

		float slope_y = deltay / (float)deltax;
		float y = (float)y0;

		for (int x = x0; x <= x1; ++x, y += slope_y)
		{
			if (x >= (int)m_width || x < 0)
				continue;
			if (y >= (int)m_height || y < 0)
				continue;

			float interp = clamp<float>((x - x0) / (float)deltax, 0, 1);

			float depth = z0 + deltaz * interp;
#if (PERSPECTIVE_CORRECTION == 1)
			float inv_w_cur = interp * delta_inv_w + inv_w0;
			depth /= inv_w_cur;
#endif

			if (depth > 1.0f || depth < 0.0f)
				continue;

			uint zVal = uint(depth * 0xFFFFFFFFUL);

			uint idx = x + (int)y * m_width;
			if (m_zBuffer[idx] < zVal)
				continue;

			Color color(*c0);
			color.R += m_ftoi(deltaC[0] * interp);
			color.G += m_ftoi(deltaC[1] * interp);
			color.B += m_ftoi(deltaC[2] * interp);
			color.A += m_ftoi(deltaC[3] * interp);

			m_bits[idx] = color.val;

			// Inversed Z in buffer
			m_zBuffer[idx] = zVal;
		}
	}
	else
	{
		if (y0 > y1)
		{
			temp = x0;
			x0 = x1;
			x1 = temp;

			temp = y0;
			y0 = y1;
			y1 = temp;

			float tempF = z0;
			z0 = z1;
			z1 = tempF;

			deltax = -deltax;
			deltay = -deltay;
			deltaz = -deltaz;

#if (PERSPECTIVE_CORRECTION == 1)
			inv_w0 = v1.inv_w;
			delta_inv_w = -delta_inv_w;
#endif

			deltaC[0] = v0.color.R - v1.color.R;
			deltaC[1] = v0.color.G - v1.color.G;
			deltaC[2] = v0.color.B - v1.color.B;
			deltaC[3] = v0.color.A - v1.color.A;
			c0 = &v1.color;
		}
		else
		{
			deltaC[0] = v1.color.R - v0.color.R;
			deltaC[1] = v1.color.G - v0.color.G;
			deltaC[2] = v1.color.B - v0.color.B;
			deltaC[3] = v1.color.A - v0.color.A;
			c0 = &v0.color;
		}

		float slope_x = deltax / (float)deltay;
		float x = (float)x0;

		for (int y = y0; y <= y1; ++y, x += slope_x)
		{
			if (x >= (int)m_width || x < 0)
				continue;
			if (y >= (int)m_height || y < 0)
				continue;

			float interp = clamp<float>((y - y0) / (float)deltay, 0, 1);

			float depth = z0 + deltaz * interp;
#if (PERSPECTIVE_CORRECTION == 1)
			float inv_w_cur = interp * delta_inv_w + inv_w0;
			depth /= inv_w_cur;
#endif

			if (depth > 1.0f || depth < 0.0f)
				continue;

			uint zVal = uint(depth * 0xFFFFFFFFUL);

			uint idx = (int)x + y * m_width;
			if (m_zBuffer[idx] < zVal)
				continue;

			Color color(*c0);
			color.R += m_ftoi(deltaC[0] * interp);
			color.G += m_ftoi(deltaC[1] * interp);
			color.B += m_ftoi(deltaC[2] * interp);
			color.A += m_ftoi(deltaC[3] * interp);

			m_bits[idx] = color.val;

			// Inversed Z in buffer
			m_zBuffer[idx] = zVal;
		}
	}
}

// See header file for the variable naming scheme
void Drawer::interpolateAndPlotScanline(
	int y,

	const Vec3 &spR,
	const Vec3 &dR1,
	const Vec3 &dRO,

	const float *inv_wR,
	float inv_wR1,
	float inv_wRO,

	const Vec4 *cR,
	const Vec4 &cR1,
	const Vec4 &cRO,

	const Vec3 *nR,
	const Vec3 &nR1,
	const Vec3 &nRO,

	const Vec3 *tcR,
	const Vec3 &tcR1,
	const Vec3 &tcRO
	)
{
	float t1 = (y - spR.y) / dR1.y;

	if (t1 < 0.0f || t1 > 1.0f)
		return;

	float x_from = t1 * dR1.x + spR.x;
	float z_from = t1 * dR1.z + spR.z;

	float t2 = (y - spR.y) / dRO.y;

	if (t2 < 0.0f || t2 > 1.0f)
		return;

	float x_to = t2 * dRO.x + spR.x;
	float z_to = t2 * dRO.z + spR.z;

	if (x_to == x_from)
		return;

	Vec3 sp_beg, sp_end, sp_diff, sp_cur;
	Vec3 *sp_from, *sp_to;

	Vec4 c_beg, c_end, c_diff, c_cur;
	Vec4 *c_from, *c_to;

	Vec3 n_beg, n_end, n_diff, n_cur;
	Vec3 *n_from, *n_to;

	Vec3 tc_beg, tc_end, tc_diff, tc_cur;
	Vec3 *tc_from, *tc_to;

#if (PERSPECTIVE_CORRECTION == 1)
	float inv_w_beg, inv_w_end, inv_w_diff, inv_w_cur;
	float *inv_w_from, *inv_w_to;

	// inv_w interpolation
	inv_w_beg = inv_wR1 * t1 + *inv_wR;
	inv_w_end = inv_wRO * t2 + *inv_wR;

	inv_w_from = &inv_w_beg;
	inv_w_to = &inv_w_end;
#endif

	// Vertex coord interpolation
	sp_beg = t1 * dR1 + spR;
	sp_end = t2 * dRO + spR;

	sp_from = &sp_beg;
	sp_to = &sp_end;

	// Color interpolation
	c_beg = t1 * cR1 + *cR;
	c_end = t2 * cRO + *cR;

	c_from = &c_beg;
	c_to = &c_end;

	// Normal interpolation
	n_beg = nR1 * t1 + *nR;
	n_end = nRO * t2 + *nR;

	n_from = &n_beg;
	n_to = &n_end;

	// Texture coordinates interpolation
	tc_beg = tcR1 * t1 + *tcR;
	tc_end = tcRO * t2 + *tcR;

	tc_from = &tc_beg;
	tc_to = &tc_end;

	if (x_to < x_from)
	{
		float temp = x_to;
		x_to = x_from;
		x_from = temp;

		temp = z_to;
		z_to = z_from;
		z_from = temp;

#if (PERSPECTIVE_CORRECTION == 1)
		inv_w_from = &inv_w_end;
		inv_w_to = &inv_w_beg;
#endif

		sp_from = &sp_end;
		sp_to = &sp_beg;

		c_from = &c_end;
		c_to = &c_beg;

		n_from = &n_end;
		n_to = &n_beg;

		tc_from = &tc_end;
		tc_to = &tc_beg;
	}

#if (PERSPECTIVE_CORRECTION == 1)
	inv_w_diff = *inv_w_to - *inv_w_from;
#endif

	sp_diff = *sp_to - *sp_from;

	c_diff = *c_to - *c_from;

	n_diff = *n_to - *n_from;

	tc_diff = *tc_to - *tc_from;

	float inv_x_diff = 1.0f / (x_to - x_from);
	float z_diff = z_to - z_from;
	for (int x = (int)x_from; x <= (int)x_to; ++x)
	{
		if (x < 0 || x >= (int)m_width)
			continue;

#if (DBG_RASTERIZE_EDGES == 1)
		if (x != (int)x_from && x != (int)x_to)
			continue;
#endif

		float interpolant = clamp<float>((x - x_from) * inv_x_diff, 0, 1);
		float depth = interpolant * z_diff + z_from;
#if (PERSPECTIVE_CORRECTION == 1)
		inv_w_cur = interpolant * inv_w_diff + *inv_w_from;
		depth /= inv_w_cur;
#endif

		if (depth < 0.0f || depth > 1.0f)
			continue;

		uint z = static_cast<uint>(depth * 0xFFFFFFFFUL);
		uint idx = x + y * m_width;

		if (m_zBuffer[idx] < z)
			continue;

		c_cur = interpolant * c_diff + *c_from;

		Vec3 n_temp = interpolant * n_diff + *n_from;

#if (PERSPECTIVE_CORRECTION == 1)
		c_cur /= inv_w_cur;
		n_temp /= inv_w_cur;
#endif

#if (NORMALIZATION == 0)
		n_cur = n_temp;
#elif (NORMALIZATION == 1)
		n_cur = n_temp.fastNormalize();
#elif (NORMALIZATION == 2)
		n_cur = n_temp;
		n_cur.normalize();
#endif

#ifdef NORMALIZATION_EPS
		float sqLen = n_cur.sqLen();
		if (fast_abs(sqLen - 1.0f) > NORMALIZATION_EPS)
		{
			__asm
			{
				int 3
			};
		}
#endif

		if (!useSphericalEnvMapping)
		{
			tc_cur = interpolant * tc_diff + *tc_from;
#if (PERSPECTIVE_CORRECTION == 1)
			tc_cur /= inv_w_cur;
#endif
		}
		else
		{
			sp_cur = interpolant * sp_diff + *sp_from;

			// Screen-space vertex coordinates
			Vec3 sp_ss = sp_cur;
			sp_ss.x -= m_width * 0.5f;
			sp_ss.y -= m_height * 0.5f;
			sp_ss.x /= m_width * 0.5f;
			sp_ss.y /= m_height * 0.5f;

			Vec3 eyePoint = Vec3C(0.0f, 0.0f, -1.0f);
			Vec3 eyeDir = (sp_ss - eyePoint).normalize();
			Vec3 r = eyeDir - 2.0f * n_cur.dot(eyeDir) * n_cur;
			r.z += 1.0f;
			float m = 1.0f / 2.0f * sqrtf(r.x*r.x + r.y*r.y + r.z*r.z) + 1e-5f;
			tc_cur = Vec3C(r.x * m + 0.5f, r.y * m + 0.5f, 0.0f);
		}

		Color out_color;
		if (tex2D)
		{
			Color texel = tex2D->get(tc_cur.x, tc_cur.y);
#if 1
			out_color.R = (uchar)(c_cur.x * texel.R);
			out_color.G = (uchar)(c_cur.y * texel.G);
			out_color.B = (uchar)(c_cur.z * texel.B);
			out_color.A = (uchar)(c_cur.w * texel.A);
#else
			out_color.R = (uchar)clamp(c_cur.x * texel.R, 0.0f, 255.0f);
			out_color.G = (uchar)clamp(c_cur.y * texel.G, 0.0f, 255.0f);
			out_color.B = (uchar)clamp(c_cur.z * texel.B, 0.0f, 255.0f);
			out_color.A = (uchar)clamp(c_cur.w * texel.A, 0.0f, 255.0f);
#endif
		}

#if (DBG_TRIANGLE_RASTER == 1)
		out_color.R = m_ftoi(n_cur.x * 127 + 128);
		out_color.G = m_ftoi(n_cur.y * 127 + 128);
		out_color.B = m_ftoi(n_cur.z * 127 + 128);
#elif (DBG_TRIANGLE_RASTER == 2)
		out_color.R = m_ftoi(tc_cur.x * 255);
		out_color.G = m_ftoi(tc_cur.y * 255);
		out_color.B = m_ftoi(tc_cur.z * 255);
#endif

		m_bits[idx] = out_color.val;
		m_zBuffer[idx] = z;
	}
}

void Drawer::drawTriangle(const Vertex &v0, const Vertex &v1, const Vertex &v2)
{
	float cz =	(v1.coord.x - v0.coord.x)*(v2.coord.y - v0.coord.y) -
				(v1.coord.y - v0.coord.y)*(v2.coord.x - v0.coord.x);

#if (BACKFACE_CULL_CCW == 1)
	// CW is rendered
	if (cz > 0.0f)
		return;
#endif

#if (BACKFACE_CULL_CW == 1)
	// CCW is rendered
	if (cz < 0.0f)
		return;
#endif

	// Sorted points
	Vec3 sp0(v0.coord);	// Top
	Vec3 sp1(v1.coord);	// Middle
	Vec3 sp2(v2.coord);	// Bottom

	// "inv_w" is not used for non-perspective-corrected rendering down the pipe
	const float *inv_w0 = &v0.inv_w, *inv_w1 = &v1.inv_w, *inv_w2 = &v2.inv_w;
	const Vec4 *c0 = &v0.fColor, *c1 = &v1.fColor, *c2 = &v2.fColor;
	const Vec3 *n0 = &v0.normal, *n1 = &v1.normal, *n2 = &v2.normal;
	const Vec3 *tc0 = &v0.texCoord, *tc1 = &v1.texCoord, *tc2 = &v2.texCoord;

#define EXCH_VERTEX_DATA(idx0, idx1) \
		sp##idx0 = (v##idx1).coord; \
		inv_w##idx0 = &(v##idx1).inv_w; \
		c##idx0 = &(v##idx1).fColor; \
		n##idx0 = &(v##idx1).normal; \
		tc##idx0 = &(v##idx1).texCoord;

	// Sort by Y
	if (v0.coord.y < v1.coord.y)
	{
		if (v0.coord.y < v2.coord.y)
		{
			// sp0 already p0
			//sp0 = p0;

			if (v1.coord.y < v2.coord.y)
			{
				// points already set this way
				//sp1 = p1;
				//sp2 = p2;
			}
			else
			{
				EXCH_VERTEX_DATA(1, 2);
				EXCH_VERTEX_DATA(2, 1);
			}
		}
		else
		{
			EXCH_VERTEX_DATA(0, 2);
			EXCH_VERTEX_DATA(1, 0);
			EXCH_VERTEX_DATA(2, 1);
		}
	}
	else
	{
		// p1 < p0
		if (v1.coord.y < v2.coord.y)
		{
			EXCH_VERTEX_DATA(0, 1);

			if (v0.coord.y < v2.coord.y)
			{
				EXCH_VERTEX_DATA(1, 0);

				// sp2 already p2
				//sp2 = p2;
			}
			else
			{
				EXCH_VERTEX_DATA(1, 2);
				EXCH_VERTEX_DATA(2, 0);
			}
		}
		else
		{
			EXCH_VERTEX_DATA(0, 2);
			// sp1 already p1
			//sp1 = p1;
			EXCH_VERTEX_DATA(2, 0);
		}
	}

#undef EXCH_VERTEX_DATA

	// Triangle is actually line
	if (int(sp0.y - sp2.y) == 0)
		return;

	Vec3 d01 = sp1 - sp0;
	Vec3 d02 = sp2 - sp0;
	Vec3 d21 = sp1 - sp2;

	Vec4 c01 = *c1 - *c0;
	Vec4 c02 = *c2 - *c0;
	Vec4 c21 = *c1 - *c2;

	// Not used if perspective correction is disabled
	float inv_w01 = *inv_w1 - *inv_w0;
	float inv_w02 = *inv_w2 - *inv_w0;
	float inv_w21 = *inv_w1 - *inv_w2;

	Vec3 n01 = *n1 - *n0;
	Vec3 n02 = *n2 - *n0;
	Vec3 n21 = *n1 - *n2;

	Vec3 tc01 = *tc1 - *tc0;
	Vec3 tc02 = *tc2 - *tc0;
	Vec3 tc21 = *tc1 - *tc2;

	int d21i = m_ftoi(d21.y);

	for (int y = (int)sp0.y; y <= (int)sp2.y; ++y)
	{
		if (y < 0 || y >= (int)m_height)
			continue;

		if (y < sp1.y || (y == sp1.y && (d21i == 0)))
		{
			interpolateAndPlotScanline(
				y,
				sp0,
				d01,
				d02,
				inv_w0,
				inv_w01,
				inv_w02,
				c0,
				c01,
				c02,
				n0,
				n01,
				n02,
				tc0,
				tc01,
				tc02
				);
		}
		else
		{
			interpolateAndPlotScanline(
				y,
				sp2,
				d21,
				-d02,
				inv_w2,
				inv_w21,
				-inv_w02,
				c2,
				c21,
				-c02,
				n2,
				n21,
				-n02,
				tc2,
				tc21,
				-tc02
				);
		}
	}
}
