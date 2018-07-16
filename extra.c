
#include "compaxx.h"
#include "compaxx_int.h"

#include <math.h>
#include <stdio.h>

float matrixDet(const Matrix* m) {
  return
    m->a1 * m->b2 * m->c3 -
    m->a1 * m->b3 * m->c2 -
    m->a2 * m->b1 * m->c3 +
    m->a2 * m->b3 * m->c1 +
    m->a3 * m->b1 * m->c2 -
    m->a3 * m->b2 * m->c1;
}

void matrixInv(const Matrix* m, Matrix* res) {
  float determ = matrixDet(m);

  float a11 = m->a1;
  float a12 = m->a2;
  float a13 = m->a3;
  float a21 = m->b1;
  float a22 = m->b2;
  float a23 = m->b3;
  float a31 = m->c1;
  float a32 = m->c2;
  float a33 = m->c3;

  res->a1 = (a22 * a33 - a23 * a32) / determ;
  res->a2 = (a13 * a32 - a12 * a33) / determ;
  res->a3 = (a12 * a23 - a13 * a22) / determ;
  res->b1 = (a23 * a31 - a21 * a33) / determ;
  res->b2 = (a11 * a33 - a13 * a31) / determ;
  res->b3 = (a13 * a21 - a11 * a23) / determ;
  res->c1 = (a21 * a32 - a22 * a31) / determ;
  res->c2 = (a12 * a31 - a11 * a32) / determ;
  res->c3 = (a11 * a22 - a12 * a21) / determ;
}

void planeFromThreePoints(const Point* p1, const Point* p2, const Point* p3, Point* cartesian) {
  Point normal;

  Point ab = {
    (float)(p2->x - p1->x),
    (float)(p2->y - p1->y),
    (float)(p2->z - p1->z)
  };

  Point ac = {
    (float)(p3->x - p1->x),
    (float)(p3->y - p1->y),
    (float)(p3->z - p1->z)
  };
  crossProduct(&ab, &ac, &normal);
  
  //normal.x = ((float)(p2->y) - (float)(p1->y)) * ((float)(p3->z) - (float)(p1->z)) - ((float)(p2->z) - (float)(p1->z)) * ((float)(p3->y) - (float)(p1->y));
  //normal.y = ((float)(p2->x) - (float)(p1->x)) * ((float)(p3->z) - (float)(p1->z)) - ((float)(p2->z) - (float)(p1->z)) * ((float)(p3->x) - (float)(p1->x));
  //normal.z = ((float)(p2->x) - (float)(p1->x)) * ((float)(p3->y) - (float)(p1->y)) - ((float)(p2->y) - (float)(p1->y)) * ((float)(p3->x) - (float)(p1->x));

  if (normal.x < 0) {
    normal.x = -normal.x;
    normal.y = -normal.y;
    normal.z = -normal.z;
  }
  //printPt(&normal, "Normal");
  
  Point p = { p1->x, p1->y, p1->z };
  normalToCartesian(&normal, &p, cartesian);
}

float rectangleArea(const Point* pt1, const Point* pt2, const Point* pt3) {
  Point ab = { pt1->x - pt2->x, pt1->y - pt2->y, pt1->z - pt2->z };
  Point ac = { pt1->x - pt3->x, pt1->y - pt3->y, pt1->z - pt3->z };
  Point cp;
  crossProduct(&ab, &ac, &cp);
  return 0.5 * (sqrt(cp.x * cp.x + cp.y * cp.y + cp.z * cp.z));
}

short finalizeCalibrationTrian(const CalibrationContext* ctx, Calibration* cal) {
  int ofs1 = 0;
  int ofs2 = ctx->pointCount / 3;
  int ofs3 = ctx->pointCount / 3 * 2;

  if (ofs1 == ofs2 || ofs2 == ofs3 || ofs3 == ofs1)
    return E_NOT_ENOUGH_CALIBRATION_POINTS;

  int idx1 = ofs1;
  int idx2 = ofs2;
  int idx3 = ofs3;

  Point accum = { 0, 0, 0 };
  float planeCount = 0;
  
  while (idx1 < ofs2 && idx2 < ofs3 && idx3 < ctx->pointCount) {
    Point plane;
    planeFromThreePoints(&(ctx->points[idx1].sensorData),
			 &(ctx->points[idx2].sensorData),
			 &(ctx->points[idx3].sensorData),
			 &plane);
    //printPt(&plane, "Current plane");

    float weight = rectangleArea(&(ctx->points[idx1].sensorData),
				 &(ctx->points[idx2].sensorData),
				 &(ctx->points[idx3].sensorData));
    accum.x += plane.x * weight;
    accum.y += plane.y * weight;
    accum.z += plane.z * weight;
    planeCount += weight;

    idx1++;
    idx2++;
    idx3++;
  }
  cal->planeA = accum.x / planeCount;
  cal->planeB = accum.y / planeCount;
  cal->planeC = accum.z / planeCount;
  return E_SUCCESS;
}
