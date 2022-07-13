#ifndef Walk_H
#define Walk_H

#include "Tools.h"
#include "LocoState.h"

class Walk
{
public:
    Walk(float mt, float mc, float, float l);
    Walk();
    ~Walk();

    LocoState vr_op(std::string op, LocoState current_state, int theta_diff, int d, int c);
    int get_dist();
    int get_cost();
    std::vector<std::string> get_op();

private:

    LocoState _null(LocoState current_state, int theta_diff);
    PhyState _null(PhyState current_state, int theta_diff);
    LocoState _reset(LocoState current_state, int theta_diff, int theta_reset);
    LocoState _translation(LocoState current_state, int theta_diff);
    LocoState _rotation(LocoState current_state, int theta_diff);
    LocoState _curvature(LocoState current_state, int theta_diff);
    unsigned _get_x_move(unsigned angle) const;
    unsigned _get_y_move(unsigned angle) const;
    LocoState _curvature_minus_45(LocoState current_state, int theta_diff);
    LocoState _curvature_45(LocoState current_state, int theta_diff);

    float mt;
    float mc;
    float mr;
    float l;
    int distance;
    int cost;
    int num_directions;

    std::map<int, int> x_move;
    std::map<int, int> y_move;
    std::vector<std::string> _vr_op_set;
    int angles[8];

    //ShaoHeng implements
    std::map<std::string, float> _vr_op_costs;
};

#endif // Walk_H
