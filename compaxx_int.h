
#ifndef __COMPAXX_INT_H__
#define __COMPAXX_INT_H__

#include "compaxx.h"

typedef struct {
  float x;
  float y;
  float z;
} Point;

typedef struct {
  float xx;
  float xy;
  float xz;
  float yy;
  float yz;
  float zz;
} CovarianceMatrix;

typedef struct {
  float a1, a2, a3, b1, b2, b3, c1, c2, c3;
} Matrix;

#ifdef NULL
#undef NULL
#endif

#define NULL (void*)0

void centroid(const CalibrationCtxPoint* points, short numPoints, Point* result);

void covariance(const CalibrationCtxPoint* points, short numPoints, const Point* centroid, CovarianceMatrix* result);

float matrixDet(const Matrix* m);

void matrixInv(const Matrix* m, Matrix* res);

void printPt(const Point* pt, const char* msg);

void planeFromThreePoints(const SensorData* p1, const SensorData* p2, const SensorData* p3, Point* cartesian);

float ptPlaneDistance(const SensorData* pt, const Calibration* plane);

void normalToCartesian(const Point* normal, const Point* pt, Point* cartesian);

#endif
