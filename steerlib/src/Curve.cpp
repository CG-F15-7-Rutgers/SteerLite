//
// Copyright (c) 2015 Mahyar Khayatkhoei
// Copyright (c) 2009-2014 Shawn Singh, Glen Berseth, Mubbasir Kapadia, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
//

#include <algorithm>
#include <vector>
#include <util/Geometry.h>
#include <util/Curve.h>
#include <util/Color.h>
#include <util/DrawLib.h>
#include "Globals.h"
#include <math.h>

using namespace Util;

Curve::Curve(const CurvePoint& startPoint, int curveType) : type(curveType)
{
	controlPoints.push_back(startPoint);
}

Curve::Curve(const std::vector<CurvePoint>& inputPoints, int curveType) : type(curveType)
{
	controlPoints = inputPoints;
	sortControlPoints();
}

// Add one control point to the vector controlPoints
void Curve::addControlPoint(const CurvePoint& inputPoint)
{
	controlPoints.push_back(inputPoint);
	sortControlPoints();
}

// Add a vector of control points to the vector controlPoints
void Curve::addControlPoints(const std::vector<CurvePoint>& inputPoints)
{
	for (int i = 0; i < inputPoints.size(); i++)
		controlPoints.push_back(inputPoints[i]);
	sortControlPoints();
}

// Draw the curve shape on screen, usign window as step size (bigger window: less accurate shape)
void Curve::drawCurve(Color curveColor, float curveThickness, int window)
{
#ifdef ENABLE_GUI
    Point startPoint;
    Point endPoint;

    // Robustness: make sure there is at least two control point: start and end points
    if (!checkRobust()) {return;}

    // Move on the curve from t=0 to t=finalPoint, using window as step size, and linearly interpolate the curve points
    if (!calculatePoint(startPoint, 0)) {return;}

    for (int t = 1; t < controlPoints.back().time; t += window)
    {
        if (!calculatePoint(endPoint, t)) {return;}
        DrawLib::drawLine(startPoint, endPoint, curveColor, curveThickness);
        startPoint = endPoint;
    }

	return;
#endif
}

bool timeSort(CurvePoint point1, CurvePoint point2)
{
    return point1.time < point2.time;
}

// Sort controlPoints vector in ascending order: min-first
void Curve::sortControlPoints()
{
    std::sort(controlPoints.begin(),controlPoints.end(), timeSort);
}

// Calculate the position on curve corresponding to the given time, outputPoint is the resulting position
bool Curve::calculatePoint(Point& outputPoint, float time)
{
	// Robustness: make sure there is at least two control point: start and end points
	if (!checkRobust())
		return false;

	// Define temporary parameters for calculation
	unsigned int nextPoint;
	float normalTime, intervalTime;

	// Find the current interval in time, supposing that controlPoints is sorted (sorting is done whenever control points are added)
	if (!findTimeInterval(nextPoint, time))
		return false;

	// Calculate position at t = time on curve
	if (type == hermiteCurve)
	{
		outputPoint = useHermiteCurve(nextPoint, time);
	}
	else if (type == catmullCurve)
	{
		outputPoint = useCatmullCurve(nextPoint, time);
	}

	// Return
	return true;
}

// Check Roboustness
bool Curve::checkRobust()
{
    if (controlPoints.size() > 2) {return true;}
	return false;
}

// Find the current time interval (i.e. index of the next control point to follow according to current time)
bool Curve::findTimeInterval(unsigned int& nextPoint, float time)
{
	for (int i = 0; i < controlPoints.size(); i++)
	{
        if (controlPoints[i].time > time)
        {
            nextPoint = i;
            return true;
        }
	}

	return false;
}

// Implement Hermite curve
Point Curve::useHermiteCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	float t, intervalTime;

	// Calculate time interval, and normal time required for later curve calculations
	intervalTime = controlPoints[nextPoint].time - controlPoints[nextPoint - 1].time;
	t = (time - controlPoints[nextPoint - 1].time) / intervalTime;

	// Calculate position at t = time on Hermite curve
	Point p0 = (2 * pow(t,3) - 3 * pow(t,2) + 1) * controlPoints[nextPoint-1].position;
	Vector t0 = (pow(t,3) - 2 * pow(t,2) + t) * intervalTime * controlPoints[nextPoint-1].tangent;
	Point p1 = (-2 * pow(t,3) + 3 * pow(t,2)) * controlPoints[nextPoint].position;
	Vector t1 = (pow(t,3) - pow(t,2)) * intervalTime * controlPoints[nextPoint].tangent;

	newPosition = p0 + t0 + p1 + t1;

	// Return result
	return newPosition;
}

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	float t, intervalTime;
	
	// Calculate time interval, and normal time required for later curve calculations
	intervalTime = controlPoints[nextPoint].time - controlPoints[nextPoint - 1].time;
	t = (time - controlPoints[nextPoint - 1].time) / intervalTime;

	// Calculate position at t = time on Catmull-Rom curve
	Point p0 = (2 * pow(t,3) - 3 * pow(t,2) + 1) * controlPoints[nextPoint-1].position;
	Point p1 = (-2 * pow(t,3) + 3 * pow(t,2)) * controlPoints[nextPoint].position;

	float b1, b2;
	Vector a1, a2;
	unsigned int nextPoint2 = nextPoint - 1;

	// First Point
	if (nextPoint2 == 0)
	{
		b1 = (controlPoints[2].time - controlPoints[0].time) / (controlPoints[2].time - controlPoints[1].time);

		a1 = (controlPoints[1].position - controlPoints[0].position) / (controlPoints[1].time - controlPoints[0].time);

		b2 = (controlPoints[0].time - controlPoints[1].time) / (controlPoints[2].time - controlPoints[1].time);

		a2 = (controlPoints[2].position - controlPoints[0].position) / (controlPoints[2].time - controlPoints[0].time);
	}
	// Last Point
	else if (nextPoint2 == controlPoints.size() - 1)
	{
		b1 = (controlPoints[nextPoint2].time-controlPoints[nextPoint2-2].time) / (controlPoints[nextPoint2-1].time-controlPoints[nextPoint2-2].time);

		a1 = (controlPoints[nextPoint2].position - controlPoints[nextPoint2-1].position) / (controlPoints[nextPoint2].time - controlPoints[nextPoint2-1].time);

		b2 = (controlPoints[nextPoint2-1].time-controlPoints[nextPoint2].time) / (controlPoints[nextPoint2-1].time-controlPoints[nextPoint2-2].time);

		a2 = (controlPoints[nextPoint2].position - controlPoints[nextPoint2-2].position) / (controlPoints[nextPoint2].time - controlPoints[nextPoint2-2].time);
	}
	// Any Other Point
	else
	{
		b1 = (controlPoints[nextPoint2].time - controlPoints[nextPoint2-1].time) / (controlPoints[nextPoint2+1].time - controlPoints[nextPoint2-1].time);

		a1 = (controlPoints[nextPoint2+1].position - controlPoints[nextPoint2].position) / (controlPoints[nextPoint2+1].time - controlPoints[nextPoint2].time);

		b2 = (controlPoints[nextPoint2+1].time - controlPoints[nextPoint2].time) / (controlPoints[nextPoint2+1].time - controlPoints[nextPoint2-1].time);

		a2 = (controlPoints[nextPoint2].position - controlPoints[nextPoint2-1].position) / (controlPoints[nextPoint2].time - controlPoints[nextPoint2-1].time);
	}

	Vector t0 = (pow(t,3) - 2 * pow(t,2) + t) * intervalTime * (b1 * a1 + b2 * a2);


	// First Point
	if (nextPoint == 0)
	{
		b1 = (controlPoints[2].time - controlPoints[0].time) / (controlPoints[2].time - controlPoints[1].time);

		a1 = (controlPoints[1].position - controlPoints[0].position) / (controlPoints[1].time - controlPoints[0].time);

		b2 = (controlPoints[0].time - controlPoints[1].time) / (controlPoints[2].time - controlPoints[1].time);

		a2 = (controlPoints[2].position - controlPoints[0].position) / (controlPoints[2].time - controlPoints[0].time);
	}
	// Last Point
	else if (nextPoint == controlPoints.size() - 1)
	{
		b1 = (controlPoints[nextPoint].time-controlPoints[nextPoint-2].time) / (controlPoints[nextPoint-1].time-controlPoints[nextPoint-2].time);

		a1 = (controlPoints[nextPoint].position - controlPoints[nextPoint-1].position) / (controlPoints[nextPoint].time - controlPoints[nextPoint-1].time);

		b2 = (controlPoints[nextPoint-1].time-controlPoints[nextPoint].time) / (controlPoints[nextPoint-1].time-controlPoints[nextPoint-2].time);

		a2 = (controlPoints[nextPoint].position - controlPoints[nextPoint-2].position) / (controlPoints[nextPoint].time - controlPoints[nextPoint-2].time);
	}
	// Any Other Point
	else
	{
		b1 = (controlPoints[nextPoint].time - controlPoints[nextPoint-1].time) / (controlPoints[nextPoint+1].time - controlPoints[nextPoint-1].time);

		a1 = (controlPoints[nextPoint+1].position - controlPoints[nextPoint].position) / (controlPoints[nextPoint+1].time - controlPoints[nextPoint].time);

		b2 = (controlPoints[nextPoint+1].time - controlPoints[nextPoint].time) / (controlPoints[nextPoint+1].time - controlPoints[nextPoint-1].time);

		a2 = (controlPoints[nextPoint].position - controlPoints[nextPoint-1].position) / (controlPoints[nextPoint].time - controlPoints[nextPoint-1].time);
	}

	Vector t1 = (pow(t,3) - pow(t,2)) * intervalTime * (b1 * a1 + b2 * a2);

	newPosition = p0 + t0 + p1 + t1;


	// Return result
	return newPosition;
}
