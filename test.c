
#include <stdio.h>
#include <assert.h>
#include "compaxx.h"
#include "compaxx_int.h"

#include <math.h>
#include <stdlib.h>

#define ASSERT_EQ(a, b, err) assert(fabs(a - b) < err)

#define RUNTEST(name) \
  if (rc == E_SUCCESS) { \
    printf(#name "\n"); \
    rc = name(); \
  }

int testCalibration(const SensorData* points, short nPoints, const Calibration* expected, float tolerance) {
  int rc = E_SUCCESS;
  CalibrationContext ctx;

  rc = startCalibration(&ctx);
  assert(rc == E_SUCCESS);
  
  int i;
  for (i=0; i<nPoints; i++) {
    rc = addCalibrationPoint(&ctx, &(points[i]), NULL);
    assert(rc == E_SUCCESS);
  }
  Calibration cal;
  float quality;
  rc = finalizeCalibration(&ctx, &cal, &quality);
  assert(rc == E_SUCCESS);

  /*
  printf("planeA (x) = %f\n", cal.planeA);
  printf("planeB (y) = %f\n", cal.planeB);
  printf("planeC (z) = %f\n", cal.planeC);
  */

  float maxErr = fmax(fmax(fabs(cal.planeA - expected->planeA), fabs(cal.planeB - expected->planeB)), fabs(cal.planeC - expected->planeC));
  //printf("Tolerance: %f\n", (float)tolerance);
  ASSERT_EQ(cal.planeA, expected->planeA, tolerance);
  ASSERT_EQ(cal.planeB, expected->planeB, tolerance);
  ASSERT_EQ(cal.planeC, expected->planeC, tolerance);

  // Now test RMSE
  float mse = 0.0;
  for (i=0; i<nPoints; i++) {
    float dist = ptPlaneDistance(&(points[i]), &cal);
    mse += dist * dist;
  }
  float rmse = sqrt(mse / (float)nPoints);
  printf("%f %f %f\n", rmse, maxErr, quality);

  return 0;
}

int testCoarseCalibrationXYplane() {
  SensorData points[] = {
    { 1, 0, 0 },
    { 0, 1, 0 },
    { 0, 0, 0 },
  };
  Calibration expected = {
    0.0, 0.0, 1.0, 0.0
  };
  return testCalibration(points, sizeof(points) / sizeof(SensorData), &expected, 0.0001);
}

float randFloat(float from, float to) {
  return (float)(rand()) * (to - from) / (float)RAND_MAX + from;
}

int testCoarseCalibrationRandomPlane() {
  float a = 0.1;
  float b = 0.2;
  float c = 0.3;
  float d = sqrt(1.0 - a*a - b*b - c*c);

  // Generate points on the plane
  #define NPOINTS 36
  #define RANGE 1000
  SensorData points[NPOINTS];
  int i;
  float noise[] = {
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 6.0, 10.0, 15.0, 20.0, 30.0, 50.0, 80.0, 100.0, 140.0, 170.0, 200.0, 300.0, -1.0
  };
  int j = 0;
  while (noise[j] >= 0.0) {
    for (i=0; i<NPOINTS; i++) {
      float x0 = randFloat(0, RANGE);
      float y0 = randFloat(0, RANGE);
      float z0 = (-a * x0 - b * y0 - d) / c;

      float nx = randFloat(-noise[j], noise[j]);
      //printf("nx: %f\n", (float)nx);
      float ny = randFloat(-noise[j], noise[j]);
      float nz = randFloat(-noise[j], noise[j]);
      points[i].x = round(x0 + nx);//round(x0);
      points[i].y = round(y0 + ny);//round(y0);
      points[i].z = round(z0 + nz);//round(z0);

      //printf("%f, %f, %f\n", (float)(points[i].x), (float)(points[i].y), (float)(points[i].z));
    }
    Calibration expected = {
      a, b, c, 0.0
    };
    float tolerance = (noise[j] == 0 ? 0.4 : log(noise[j]) / 20.0 + 0.4);
    printf("%f %f ", (float)(noise[j]), tolerance);
    testCalibration(points, NPOINTS, &expected, tolerance);
    j++;
  }
  return E_SUCCESS;
}

int testMatrixInv() {
  Matrix m = {
    1, 2, 3,
    0, 1, 4,
    5, 6, 0
  };

  Matrix inverted;

  matrixInv(&m, &inverted);
  printf("%f, %f, %f\n", inverted.a1, inverted.a2, inverted.a3);
  printf("%f, %f, %f\n", inverted.b1, inverted.b2, inverted.b3);
  printf("%f, %f, %f\n", inverted.c1, inverted.c2, inverted.c3);

  ASSERT_EQ(inverted.a1, -24, 0.0001);
  ASSERT_EQ(inverted.a2, 18, 0.0001);
  ASSERT_EQ(inverted.a3, 5, 0.0001);
  ASSERT_EQ(inverted.b1, 20, 0.0001);
  ASSERT_EQ(inverted.b2, -15, 0.0001);
  ASSERT_EQ(inverted.b3, -4, 0.0001);
  ASSERT_EQ(inverted.c1, -5, 0.0001);
  ASSERT_EQ(inverted.c2, 4, 0.0001);
  ASSERT_EQ(inverted.c3, 1, 0.0001);
  return E_SUCCESS;
}

int testPlaneFromThreePoints() {
  float a = 0.1;
  float b = 0.2;
  float c = 0.3;
  float d = sqrt(1.0 - a*a - b*b - c*c);

  printf("d = %f\n", d);
  // Generate points on the plane
  #define TNPOINTS 3
  #define TRANGE 1000
  SensorData points[TNPOINTS];
  int i;
  for (i=0; i<TNPOINTS; i++) {
    float x0 = (float)(rand()) * TRANGE / (float)RAND_MAX;;
    float y0 = (float)(rand()) / (float)RAND_MAX * TRANGE;
    float z0 = (-a * x0 - b * y0 - d) / c;
    points[i].x = round(x0);
    points[i].y = round(y0);
    points[i].z = round(z0);

    printf("%i, %i, %i\n", (int)(points[i].x), (int)(points[i].y), (int)(points[i].z));
  }
  Calibration expected = {
    a, b, c, 0.0
  };
  Point plane;
  planeFromThreePoints(&points[0], &points[1], &points[2], &plane);
  printPt(&plane, "Plane");
  ASSERT_EQ(plane.x, a, 0.4);
  ASSERT_EQ(plane.y, b, 0.4);
  ASSERT_EQ(plane.z, c, 0.4);
  return E_SUCCESS;
}

int main(int argc, char** argv) {
  int rc = E_SUCCESS;

  srand(123); // Fix randomness for repeatability

  RUNTEST(testCoarseCalibrationXYplane);
  RUNTEST(testCoarseCalibrationRandomPlane);
  RUNTEST(testPlaneFromThreePoints);
  RUNTEST(testMatrixInv);

  if (rc == E_SUCCESS)
    printf("SUCCESS\n");
  else
    printf("** FAILED\n");
  return 0;
}
