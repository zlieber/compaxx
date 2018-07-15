
#include "compaxx.h"
#include "compaxx_int.h"

#include <math.h>
#include <stdio.h>

void printPt(const Point* pt, const char* msg) {
  printf("%s: (%f, %f, %f)\n", msg, pt->x, pt->y, pt->z);
}

float sq(float n) {
  return n * n;
}

void centroid(const CalibrationCtxPoint* points, short numPoints, Point* result) {
  int i;
  //printf("NumPoints: %i\n", (int)numPoints);
  result->x = result->y = result->z = 0.0;
  for (i=0; i<numPoints; i++) {
    result->x += points[i].sensorData.x;
    result->y += points[i].sensorData.y;
    result->z += points[i].sensorData.z;
  }

  //printPt(result, "intermediate");
  
  result->x /= (float)numPoints;
  result->y /= (float)numPoints;
  result->z /= (float)numPoints;
}

void covariance(const CalibrationCtxPoint* points, short numPoints, const Point* centroid, CovarianceMatrix* result) {
  result->xx = result->xy = result->xz = result->yy = result->yz = result->zz = 0.0;
  int i;
  for (i=0;  i<numPoints; i++) {
    Point r = {
      (float)(points[i].sensorData.x) - centroid->x,
      (float)(points[i].sensorData.y) - centroid->y,
      (float)(points[i].sensorData.z) - centroid->z,
    };
    result->xx += r.x * r.x;
    result->xy += r.x * r.y;
    result->xz += r.x * r.z;
    result->yy += r.y * r.y;
    result->yz += r.y * r.z;
    result->zz += r.z * r.z;
  }

  result->xx /= (float)numPoints;
  result->xy /= (float)numPoints;
  result->xz /= (float)numPoints;
  result->yy /= (float)numPoints;
  result->yz /= (float)numPoints;
  result->zz /= (float)numPoints;
}

float dotProduct(const Point* a, const Point* b) {
  return a->x * b->x + a->y * b->y + a->z * b->z;
}

void addTo(Point* a, const Point* b) {
  a->x += b->x;
  a->y += b->y;
  a->z += b->z;
}

void mulByScalar(Point* a, float scalar) {
  a->x *= scalar;
  a->y *= scalar;
  a->z *= scalar;
}

float floatAbs(float f) {
  if (f >= 0.0)
    return f;
  else
    return -f;
}

float floatMax(float a, float b) {
  return (a > b ? a : b);
}

void normalize(Point* pt) {
  // Divide by largest value first to avoid overflow.
  float largest = floatMax(floatMax(floatAbs(pt->x), floatAbs(pt->y)), floatAbs(pt->z));
  pt->x /= largest;
  pt->y /= largest;
  pt->z /= largest;
  //printPt(pt, "In normalize");
  float length = sqrt(pt->x * pt->x + pt->y * pt->y + pt->z * pt->z);
  //printf("Length: %f\n", length);
  pt->x /= length;
  pt->y /= length;
  pt->z /= length;
}

short weightedDir(const CovarianceMatrix* covar, Point* weighted_dir) {
  weighted_dir->x = weighted_dir->y = weighted_dir->z = 0.0;

  float det_x = covar->yy * covar->zz - covar->yz * covar->yz;
  float det_y = covar->xx * covar->zz - covar->xz * covar->xz;
  float det_z = covar->xx * covar->yy - covar->xy * covar->xy;

  float scaling = floatMax(floatMax(floatAbs(det_x), floatAbs(det_y)), floatAbs(det_z));

  {
    Point axis_dir = {
      det_x,
      covar->xz * covar->yz - covar->xy * covar->zz,
      covar->xy * covar->yz - covar->xz * covar->yy
    };
    float weight = det_x * det_x;
    if (dotProduct(weighted_dir, &axis_dir) < 0.0)
      weight = -weight;
    mulByScalar(&axis_dir, weight);
    addTo(weighted_dir,&axis_dir);
  }

  //printPt(weighted_dir, "weighted_dir");
  
  {
    Point axis_dir = {
      covar->xz * covar->yz - covar->xy * covar->zz,
      det_y,
      covar->xy * covar->xz - covar->yz * covar->xx
    };
    float weight = det_y * det_y;
    if (dotProduct(weighted_dir, &axis_dir) < 0.0)
      weight = -weight;
    mulByScalar(&axis_dir, weight);
    addTo(weighted_dir,&axis_dir);
  }

  //printPt(weighted_dir, "weighted_dir");

  {
    Point axis_dir = {
      covar->xy * covar->yz - covar->xz * covar->yy,
      covar->xy * covar->xz - covar->yz * covar->xx,
      det_z,
    };
    float weight = det_z * det_z;
    if (dotProduct(weighted_dir, &axis_dir) < 0.0)
      weight = -weight;
    mulByScalar(&axis_dir, weight);
    addTo(weighted_dir,&axis_dir);
  }

  //printPt(weighted_dir, "weighted_dir");
  normalize(weighted_dir);

  //printPt(weighted_dir, "weighted_dir");
  return E_SUCCESS;
}

void normalToCartesian(const Point* normal, const Point* pt, Point* cartesian) {
  float a, b, c, d;

  a = normal->x;
  b = normal->y;
  c = normal->z;
  d = -(a * pt->x + b * pt->y + c * pt->z);

  //printf("D = %f\n", d);
  float length = sqrt(a*a + b*b + c*c + d*d);
  //printf("Length: %f\n", length);
  cartesian->x = a / length;
  cartesian->y = b / length;
  cartesian->z = c / length;
}

short getHeading(const Calibration* cal, const SensorData* sensorData, float* heading) {
  return E_SUCCESS;
}

short startCalibration(CalibrationContext* ctx) {
  ctx->pointCount = 0;
  return E_SUCCESS;
}

short addCalibrationPoint(CalibrationContext* ctx, const SensorData* sensorData, const float* magneticHeading) {
  if (ctx->pointCount == MAX_SENSOR_POINTS)
    return E_TOO_MANY_COARSE_POINTS;
  ctx->points[ctx->pointCount].sensorData = *sensorData;
  if (magneticHeading == NULL)
    ctx->points[ctx->pointCount].magneticHeading = -1;
  else
    ctx->points[ctx->pointCount].magneticHeading = *magneticHeading;
  ctx->pointCount++;
  return E_SUCCESS;
}

float ptPlaneDistance(const SensorData* pt, const Calibration* plane) {
  float planeD = sqrt(1.0 - plane->planeA * plane->planeA - plane->planeB * plane->planeB - plane->planeC * plane->planeC);
  return
    floatAbs(plane->planeA * pt->x + plane->planeB * pt->y + plane->planeC * pt->z) /
    sqrt(plane->planeA * plane->planeA + plane->planeB * plane->planeB + plane->planeC * plane->planeC);
}

float rmse(const CalibrationContext* ctx, const Calibration* cal) {
  float mse = 0.0;
  int i;
  for (i=0; i<ctx->pointCount; i++) {
    float dist = ptPlaneDistance(&(ctx->points[i].sensorData), cal);
    mse += dist * dist;
  }
  return sqrt(mse / (float)(ctx->pointCount));
}

float vecLength(const SensorData* sensor) {
  return sqrt(sq((float)(sensor->x)) + sq((float)(sensor->y)) + sq((float)(sensor->z)));
}

float meanLength(const CalibrationContext* ctx) {
  int i;
  float totalLength = 0.0;
  for (i=0; i<ctx->pointCount; i++) {
    totalLength += vecLength(&(ctx->points[i].sensorData));
  }
  return totalLength / ctx->pointCount;
}

short finalizeCalibration(const CalibrationContext* ctx, Calibration* cal, float* quality) {
  Point centroidPt;
  CovarianceMatrix covar;
  Point normal;
  Point cartesian;
  centroid(ctx->points, ctx->pointCount, &centroidPt);
  //printPt(&centroidPt, "centroid");
  covariance(ctx->points, ctx->pointCount, &centroidPt, &covar);

  weightedDir(&covar, &normal);
  //printPt(&normal, "Normal");
  normalToCartesian(&normal, &centroidPt, &cartesian);
  cal->planeA = cartesian.x;
  cal->planeB = cartesian.y;
  cal->planeC = cartesian.z;

  if (quality)
    *quality = 100.0 - rmse(ctx, cal) / meanLength(ctx) * 100;
  return E_SUCCESS;
}
