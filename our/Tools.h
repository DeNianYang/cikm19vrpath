#ifndef TOOL_H
#define TOOL_H

#include <vector>
#include <map>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <set>
#include <algorithm>
#include <tuple>
#include <sys/time.h>

const float TRANSLATION_2_COST = 1;
const float CURVATURE_45_COST = 1;
const float ROTATION_HALF_COST = 1;
const float RESET_180_COST = 5;

extern float RESET_COST;

const float TRANS_LB = 0.86;
const float TRANS_UB = 1.26;
const float ROTATE_LB = 0.67;
const float ROTATE_UB = 1.24;

#define arraysize(x) (sizeof(x)/sizeof(*x))

struct Point
{
    Point() : x(), y() {}
    Point(int a, int b) : x(a), y(b) {}
    int x;
    int y;
    friend bool operator<(const Point& l, const Point& r)
    {
        return std::tie(l.x, l.y) < std::tie(r.x, r.y);
    }

    friend bool operator==(const Point& l, const Point& r)
    {
        return std::tie(l.x, l.y) == std::tie(r.x, r.y);
    }
};

struct FloatPoint
{
    FloatPoint() : x(), y() {}
    FloatPoint(float a, float b) : x(a), y(b) {}
    float x;
    float y;
    friend bool operator<(const FloatPoint& l, const FloatPoint& r)
    {
        return std::tie(l.x, l.y) < std::tie(r.x, r.y);
    }

    friend bool operator==(const FloatPoint& l, const FloatPoint& r)
    {
        return std::tie(l.x, l.y) == std::tie(r.x, r.y);
    }
};

struct Step
{
    Step() : delta_x(), delta_y() {}
    Step(int a, int b) : delta_x(a), delta_y(b) {}
    int delta_x;
    int delta_y;

    friend bool operator<(const Step& l, const Step& r)
    {
        return std::tie(l.delta_x, l.delta_y) < std::tie(r.delta_x, r.delta_y);
    }

    friend bool operator==(const Step& l, const Step& r)
    {
        return std::tie(l.delta_x, l.delta_y) == std::tie(r.delta_x, r.delta_y);
    }
};

struct Edge
{
    Edge(unsigned p, unsigned l) : point(p), length(l) {}
    unsigned point;
    unsigned length;
};

typedef std::vector<unsigned> kPath;
typedef std::vector<kPath> kPathes;
typedef std::map<unsigned, kPathes> kLabel;

typedef std::vector<Point> Path;
typedef std::vector<Path> Pathes;
typedef std::map<unsigned, Pathes> Label;
typedef std::vector<Point> VirPath;
typedef std::vector<VirPath> VirPathes;

struct PhyPath {

    PhyPath() : physical_path(), vr_ops() {
        distance = std::numeric_limits<int>::max();
        cost = std::numeric_limits<int>::max();
    }

    Path physical_path;
    std::vector<std::string> vr_ops;
    int distance;
    int cost;
};


typedef std::vector<PhyPath> PhyPathes;
typedef std::map<Step, Pathes> VirLabel;
typedef std::map<Step, PhyPathes> PhyLabel;

struct PhyState
{
    PhyState() : x(), y(), theta() {}
    PhyState(int a, int b, int c) : x(a), y(b), theta(c) {}
    int x;
    int y;
    int theta;

    friend bool operator<(const PhyState& l, const PhyState& r)
    {
        return std::tie(l.x, l.y, l.theta) < std::tie(r.x, r.y, r.theta);
    }
};

struct DistancePath
{
    DistancePath() : distance(0), path() {}
    DistancePath(unsigned d, const kPath& p) : distance(d), path(p) {}
    DistancePath(unsigned p) : distance(0), path(1, p) {}
    unsigned distance;
    kPath path;
    friend bool operator<(const DistancePath& l, const DistancePath& r)
    {
        return std::tie(l.distance, l.path) > std::tie(r.distance, r.path);
    }
    friend bool operator==(const DistancePath& l, const DistancePath& r)
    {
        return std::tie(l.distance, l.path) == std::tie(r.distance, r.path);
    }
};

typedef std::vector<DistancePath> DistancePathes;

struct loco_info
{
    int gamma_v_x;
    int gamma_v_y;
    int theta_v;
    int gamma_p_x;
    int gamma_p_y;
    int theta_p;
};

int arraySize(int[]);

double GetCurrentTimeSec();

struct RW
{
    float cost;
    bool reset;
    float rotate_gain;
    float translation_gain;
    bool clockwise;
};

struct Optimized_Step
{
    float cost;
    float theta;
};

#endif