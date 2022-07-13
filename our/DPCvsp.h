#ifndef DPCvsp_H
#define DPCvsp_H

#include "Tools.h"
#include "Walk.h"
#include "LocoState.h"

//#define  DEBUG

struct VRTuple
{
    VRTuple() : s(), d() {}
    VRTuple(LocoState a, int b) : s(a), d(b) {}
    LocoState s;
    int d;
    friend bool operator==(const VRTuple& l, const VRTuple& r)
    {
        return std::tie(l.s, l.d) == std::tie(r.s, r.d);
    }

    friend bool operator<(const VRTuple& l, const VRTuple& r)
    {
        return std::tie(l.s, l.d) < std::tie(r.s, r.d);
    }
};

struct VRPath
{
    std::vector<LocoState> vr_path;
    std::vector<std::string> vr_ops;
    unsigned distance;
    unsigned cost;
};

class DPCvsp
{
public:
    DPCvsp();
    ~DPCvsp();
    DPCvsp(int vx, int vy, int th_v, int px, int py, int th_p, int gx, int gy, unsigned c);

    void read_vir_map(const char* graph_file);
    void read_phy_map(const char* graph_file);
    void print_map(std::string space);
    void DP_CVSP();
    std::pair<float, float> print_path();

private:

    void _read_phy_map(const char* graph_file);
    void _read_vir_map(const char* graph_file);
    bool _is_goal(LocoState s);
    bool _is_valid_op(int theta_diff, std::string op);
    unsigned _get_cost(unsigned dist, LocoState s);
    void _gen_vr_path(std::map<VRTuple, VRTuple>& parent_state, std::map<VRTuple, std::string>& parent_op, LocoState s_t);
    bool _is_valid(LocoState s);
    int _get_type(const Point& p, std::map<Point, int> point_map) const;

    std::map<unsigned, std::map<LocoState, unsigned>> _sigma;
    int _gamma_s_v_x;
    int _gamma_s_v_y;
    int _gamma_s_p_x;
    int _gamma_s_p_y;
    int _theta_p;
    int _theta_v;
    int _gamma_t_v_x;
    int _gamma_t_v_y;
    unsigned _cost_bound;
    int num_directions;
    int _angles[8];
    std::map<Point, int> _vir_point_map;
    unsigned _vir_length;
    unsigned _vir_width;
    unsigned _phy_length;
    unsigned _phy_width;
    std::map<Point, int> _phy_point_map;
    std::map<unsigned, std::map<LocoState, unsigned>> sigma;
    VRPath opt_path;

};

#endif
