#include "Walk.h"
#include "Tools.h"


Walk::Walk()
{
    mt = 2;
    mc = 0.875;
    mr = 0.5;
    l = 1;
    num_directions = 8;

    for (int i = 0; i < num_directions; i++) {
        angles[i] = 45 * i;
    }

    x_move = { {0, 1}, {45, 1}, {90, 0}, {135, -1},
              {180, -1}, {225, -1}, {270, 0}, {315, 1} };
    y_move = { { 0, 0 },{ 45, 1 },{ 90, 1 },{ 135, 1 },
            { 180, 0 },{ 225, -1 },{ 270, -1 },{ 315, -1 } };

    _vr_op_set.push_back("null");
    _vr_op_set.push_back("translation 2");
    _vr_op_set.push_back("reset 90");
    _vr_op_set.push_back("reset 135");
    _vr_op_set.push_back("reset 180");
    _vr_op_set.push_back("reset 225");
    _vr_op_set.push_back("reset 270");

    _vr_op_set.push_back("rotation 0.5");
    _vr_op_set.push_back("curvature -45");
    _vr_op_set.push_back("curvature 45");

    _vr_op_costs["null"] = 0;
    _vr_op_costs["translation 2"] = TRANSLATION_2_COST;
    _vr_op_costs["reset 90"] = RESET_180_COST;
    _vr_op_costs["reset 135"] = RESET_180_COST;
    _vr_op_costs["reset 180"] = RESET_180_COST;
    _vr_op_costs["reset 225"] = RESET_180_COST;
    _vr_op_costs["reset 270"] = RESET_180_COST;
    _vr_op_costs["rotation 0.5"] = ROTATION_HALF_COST;
    _vr_op_costs["curvature -45"] = CURVATURE_45_COST;
    _vr_op_costs["curvature 45"] = CURVATURE_45_COST;
}

int Walk::get_dist() {
    return distance;
}

int Walk::get_cost() {
    return cost;
}

std::vector<std::string> Walk::get_op() {
    return _vr_op_set;
}

unsigned Walk::_get_x_move(unsigned angle) const
{
    auto ite = x_move.find(angle);
    if (ite == x_move.end())
    {
        return std::numeric_limits<unsigned>::max();
    }
    return ite->second;
}

unsigned Walk::_get_y_move(unsigned angle) const
{
    auto ite = y_move.find(angle);
    if (ite == y_move.end())
    {
        return std::numeric_limits<unsigned>::max();
    }
    return ite->second;
}


LocoState Walk::vr_op(std::string op, LocoState current_state, int theta_diff, int d, int c) {
    LocoState next_state = LocoState();

    if (op == "null") {
        distance = d + 1;
        //cost = c;
        next_state = _null(current_state, theta_diff);
    }
    else if (op == "translation 2") {
        distance = d + 1 * mt;
        //cost = c + TRANSLATION_2_COST;
        next_state = _translation(current_state, theta_diff);
    }
    else if (op == "curvature -45") {
        distance = d + 1;
        //cost = c + CURVATURE_45_COST;
        next_state = _curvature_minus_45(current_state, theta_diff);
    }
    else if (op == "curvature 45") {
        distance = d + 1;
        //cost = c + CURVATURE_45_COST;
        next_state = _curvature_45(current_state, theta_diff);
    }
    else if (op == "rotation 0.5") {
        distance = d + 1;
        //cost = c + ROTATION_HALF_COST;
        next_state = _rotation(current_state, theta_diff);
    }
    else if (op == "reset 90") {
        distance = d + 1;
        //cost = c + RESET_180_COST;
        next_state = _reset(current_state, theta_diff, 90);
    }
    else if (op == "reset 135") {
        distance = d + 1;
        //cost = c + RESET_180_COST;
        next_state = _reset(current_state, theta_diff, 135);
    }
    else if (op == "reset 180") {
        distance = d + 1;
        //cost = c + RESET_180_COST;
        next_state = _reset(current_state, theta_diff, 180);
    }
    else if (op == "reset 225") {
        distance = d + 1;
        //cost = c + RESET_180_COST;
        next_state = _reset(current_state, theta_diff, 225);
    }
    else if (op == "reset 270") {
        distance = d + 1;
        //cost = c + RESET_180_COST;
        next_state = _reset(current_state, theta_diff, 270);
    }

    cost = c + _vr_op_costs[op]; // ShaoHeng implements; increments the corresponding cost

    return next_state;
};

LocoState Walk::_null(LocoState current_state, int theta_diff) {
    LocoState next_state = current_state;
    next_state.sumAngle(theta_diff, theta_diff);
    next_state.sumPosition(x_move[next_state.getThetaV()] * l, y_move[next_state.getThetaV()] * l,
                           x_move[next_state.getThetaP()] * l, y_move[next_state.getThetaP()] * l);
    return next_state;
};


LocoState Walk::_reset(LocoState current_state, int theta_diff, int theta_reset) {
    LocoState next_state = current_state;
    next_state.sumAngle(theta_diff, theta_diff + theta_reset);
    next_state.sumPosition(x_move[next_state.getThetaV()] * l, y_move[next_state.getThetaV()] * l,
                           x_move[next_state.getThetaP()] * l, y_move[next_state.getThetaP()] * l);
    return next_state;
};

LocoState Walk::_translation(LocoState current_state, int theta_diff) {
    LocoState next_state = current_state;
    next_state.sumPosition(x_move[next_state.getThetaV()] * l * mt, y_move[next_state.getThetaV()] * l * mt,
                           x_move[next_state.getThetaP()] * l, y_move[next_state.getThetaP()] * l);
    return next_state;
};


LocoState Walk::_rotation(LocoState current_state, int theta_diff) {
    LocoState next_state = current_state;
    next_state.sumAngle(theta_diff, theta_diff * mr);
    next_state.sumPosition(x_move[next_state.getThetaV()] * l, y_move[next_state.getThetaV()] * l,
                           x_move[next_state.getThetaP()] * l, y_move[next_state.getThetaP()] * l);
    return next_state;
};

LocoState Walk::_curvature_minus_45(LocoState current_state, int theta_diff) {
    LocoState next_state = current_state;
    next_state.sumAngle(0, -45);
    next_state.sumPosition(x_move[next_state.getThetaV()] * l, y_move[next_state.getThetaV()] * l,
                           x_move[next_state.getThetaP()] * l, y_move[next_state.getThetaP()] * l);
    return next_state;
};

LocoState Walk::_curvature_45(LocoState current_state, int theta_diff) {
    LocoState next_state = current_state;
    next_state.sumAngle(0, 45);
    next_state.sumPosition(x_move[next_state.getThetaV()] * l, y_move[next_state.getThetaV()] * l,
        x_move[next_state.getThetaP()] * l, y_move[next_state.getThetaP()] * l);
    return next_state;
};

Walk::~Walk()
{
}
