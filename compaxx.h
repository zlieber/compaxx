
#ifndef __COMPAXX_H__
#define __COMPAXX_H__

// The Compact Compass Library (CompaxxLib)

#define MAX_CALIBRATION_POINTS    36
#define MAX_SENSOR_POINTS         200

typedef struct {
  float compassHeading;
  float magneticHeading;
} CalibrationPoint;

typedef struct {
  float planeA;
  float planeB;
  float planeC;

  float calibrationData[MAX_CALIBRATION_POINTS];
  int pointCount;
} Calibration;

typedef struct {
  short x;
  short y;
  short z;
} SensorData;

typedef struct {
  SensorData sensorData;
  float magneticHeading;
} CalibrationCtxPoint;

typedef struct {
  CalibrationCtxPoint points[MAX_SENSOR_POINTS];
  int pointCount;
} CalibrationContext;

#define E_SUCCESS                          0
#define E_NEED_COARSE_CALIBRATION         -1
#define E_NOT_ENOUGH_CALIBRATION_POINTS   -2
#define E_TOO_MANY_COARSE_POINTS          -3
#define E_TOO_MANY_FINE_POINTS            -4

/**
 * Returns current compass or magnetic heading, given 3-axis sensor
 * reading. Compass heading is heading read off the compass without
 * accounting for effects on the instrument of various metallic
 * objects in its surroundings. As such, this function requires that a
 * coarse calibration has been performed, so that we know compass
 * orientation relative to horizontal plane. It also requires at least
 * one fine calibration point, so that it can determine which way the
 * compass is pointing on a vessel or drone.
 *
 * Magnetic heading requires more fine calibration points (up to 36,
 * but at least 8 is recommended), and corrects the compass heading
 * for effects of environment on the instrument. It is therefore more
 * precise than a compass heading.
 *
 * @param cal Existing calibration structure. At least coarse
 *        calibration must have been performed.
 * @param sensorData 3-axis sensor data provided by the instrument.
 * @param compassHeading Pointer to variable to store result in.
 * @return Error code.
 */
short getHeading(const Calibration* cal, const SensorData* sensorData, float* heading);

/**
 * Begins the process of calibrating the instrument.
 *
 * @param ctx Pointer to existing CalibrationContext structure. There
 * is no need to initialize anything in this structure, it will be
 * done by this function.
 * @return E_SUCCESS
 */
short startCalibration(CalibrationContext* ctx);

/**
 * Adds a calibration point to the calibration context.
 *
 * Calibration points can be of two types: coarse and fine calibration
 * points.
 *
 * Coarse calibration points contain instrument reading, and nothing
 * else (magneticHeading == null). These points are used to establish
 * instrument orientation relative to the horizontal plane.
 *
 * Fine calibration points contain correct magnetic heading, as read
 * by hand bearing compass that is placed away from mangetic
 * interference.
 *
 * For successful calibration process, it is required to provide at
 * least 10 coarse calibration points, and one fine calibration
 * point. This will enable reading compass heading via getHeading.
 *
 * If, in addition to that, at leat 8 more fine calibration points are
 * provided, then getHeading will start returning a more precise
 * magnetic heading instead.
 *
 * @param ctx Existing calibration context
 * @param sensorData Sensor data as returned by the 3-axis instrument
 * @param magneticHeading Optional - if this is a fine calibration
 * point, then it should contain a correct magnetic heading as read by
 * a hand bearing compass.
 * @return Error code
 */
short addCalibrationPoint(CalibrationContext* ctx, const SensorData* sensorData, const float* magneticHeading);

/**
 * Finalizes the calibration process and initializes the Calibration
 * structure to be used in future reference.
 *
 * After this function returns, the CalibrationContext is no longer
 * needed, and the scope in which it is defined can be exited.
 *
 * @param ctx Existing calibration context
 * @param cal Points to Calibration structure. There is no need to
 * initialize it.
 * @param quality Optional quality output. Will contain quality metric
 * (in percent) of the calibration performed.
 * @return Error code.
 */
short finalizeCalibration(const CalibrationContext* ctx, Calibration* cal, float* quality);

#endif
