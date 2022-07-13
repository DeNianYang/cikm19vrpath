#include "DPCvsp.h"

DPCvsp::DPCvsp()
{
    num_directions = 8;
    for (int i = 0; i < num_directions; i++) {
        _angles[i] = 45 * i;
    }
}

DPCvsp::DPCvsp(int vx, int vy, int th_v, int px, int py, int th_p, int gx, int gy, unsigned c)
{
    _gamma_s_v_x = vx;
    _gamma_s_v_y = vy;
    _gamma_s_p_x = px;
    _gamma_s_p_y = py;
    _theta_p = th_p;
    _theta_v = th_v;
    _gamma_t_v_x = gx;
    _gamma_t_v_y = gy;
    _cost_bound = c;
    num_directions = 8;
    for (int i = 0; i < num_directions; i++) {
        _angles[i] = 45 * i;
    }
}

void DPCvsp::read_vir_map(const char* graph_file)
{
    std::string s(graph_file);
    if (s.empty())
    {
        return;
    }

    _read_vir_map(graph_file);
}

void DPCvsp::read_phy_map(const char* graph_file)
{
    std::string s(graph_file);
    if (s.empty())
    {
        return;
    }

    _read_phy_map(graph_file);
}

void DPCvsp::_read_vir_map(const char* graph_file)
{
    std::ifstream infile(graph_file);
    std::string line;
    unsigned type = 3;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        int x, y;
        if (!(iss >> x >> y)) {
            if (!(iss.str() == "pois" || iss.str() == "obs"
                || iss.str() == "width" || iss.str() == "length" || iss.str() == "refs")) {
                std::cout << "_read_graph::tags error" << "\n";
                break;
            }// error
        }

        if (iss.str() == "obs") {
            type = 0;
            continue;
        }
        else if (iss.str() == "pois") {
            type = 1;
            continue;
        }
        else if (iss.str() == "refs") {
            type = 2;
            continue;
        }
        else if (iss.str() == "length") {
            std::getline(infile, line);
            std::istringstream iss(line);
            if (!(iss >> _vir_length)) {
                std::cout << "_read_graph::length error" << "\n";
            }
            continue;
        }
        else if (iss.str() == "width") {
            std::getline(infile, line);
            std::istringstream iss(line);
            if (!(iss >> _vir_width)) {
                std::cout << "_read_graph::width error" << "\n";
            }
            continue;
        }
        _vir_point_map[Point(x, y)] = type;
    }

    for (int x = 0; x < _vir_length; x++) {
        //std::cout << x << "\n";
        for (int y = 0; y < _vir_width; y++) {
            Point gamma = Point(x, y);
            if (_get_type(gamma, _vir_point_map) == 3) { //is not obs
                _vir_point_map[Point(x, y)] = 3; //3: free space
            }
        }
    }
}

void DPCvsp::_read_phy_map(const char* graph_file)
{
    std::ifstream infile(graph_file);
    std::string line;
    unsigned type = 3;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        int x, y;
        if (!(iss >> x >> y)) {
            if (!(iss.str() == "pois" || iss.str() == "obs"
                || iss.str() == "width" || iss.str() == "length" || iss.str() == "refs")) {
                std::cout << "_read_graph::tags error" << "\n";
                break;
            }// error
        }

        if (iss.str() == "obs") {
            type = 0;
            continue;
        }
        else if (iss.str() == "pois") {
            type = 1;
            continue;
        }
        else if (iss.str() == "refs") {
            type = 2;
            continue;
        }
        else if (iss.str() == "length") {
            std::getline(infile, line);
            std::istringstream iss(line);
            if (!(iss >> _phy_length)) {
                std::cout << "_read_graph::length error" << "\n";
            }
            continue;
        }
        else if (iss.str() == "width") {
            std::getline(infile, line);
            std::istringstream iss(line);
            if (!(iss >> _phy_width)) {
                std::cout << "_read_graph::width error" << "\n";
            }
            continue;
        }
        _phy_point_map[Point(x, y)] = type;
    }

    for (int x = 0; x < _phy_length; x++) {
        //std::cout << x << "\n";
        for (int y = 0; y < _phy_width; y++) {
            Point gamma = Point(x, y);
            if (_get_type(gamma, _phy_point_map) == 3) { //is not obs
                _phy_point_map[Point(x, y)] = 3; //3: free space
            }
        }
    }
}

int DPCvsp::_get_type(const Point& p, std::map<Point, int> point_map) const
{
    auto ite = point_map.find(p);
    if (ite == point_map.end())
    {
        return  3;
    }
    return ite->second;
}

void DPCvsp::print_map(std::string space)
{
    if (space == "vir") {
        for (std::map<Point, int>::iterator it = _vir_point_map.begin(); it != _vir_point_map.end(); ++it) {
            printf("(%d,%d), %d\n", it->first.x, it->first.y, it->second);
        }
    }
    else {
        for (std::map<Point, int>::iterator it = _phy_point_map.begin(); it != _phy_point_map.end(); ++it) {
            printf("(%d,%d), %d\n", it->first.x, it->first.y, it->second);
        }
    }
    
}

void DPCvsp::DP_CVSP() {
    LocoState init_state(_gamma_s_v_x, _gamma_s_v_y, _theta_v, _gamma_s_p_x, _gamma_s_p_y, _theta_p);
    if (_is_valid(init_state)) {
        Walk one_step = Walk();
        std::vector<std::string> ops = one_step.get_op();
        std::map<VRTuple, VRTuple> parent_state;
        std::map<VRTuple, std::string> parent_op;
        

        //initialization
        for (int i = 0; i < arraySize(_angles); i++) {
            LocoState cs(_gamma_s_v_x, _gamma_s_v_y, _angles[i], _gamma_s_p_x, _gamma_s_p_y, _theta_p);
            sigma[0][cs] = 0;
            parent_state[VRTuple(cs, 0)] = VRTuple(cs, -1);
            parent_op[VRTuple(cs, 0)] = "null";
        }

        bool complete_flag = false;
        LocoState final_state;
        for (int pre_d = 0; pre_d < _vir_length * _vir_width; pre_d++) {
            //std::cout << pre_d << "\n";
            std::map<LocoState, unsigned> current_states = sigma[pre_d];

            if (current_states.size() == 0) {
#ifdef DEBUG
                std::cout << "[WARNNING] no next states" << std::endl;
#endif
                break;
            }

            for (auto it = current_states.begin(); it != current_states.end(); ++it)
            {
                LocoState s = it->first;
                if (_is_goal(s)) {
                    complete_flag = true;
                    opt_path.distance = pre_d;
                    opt_path.cost = _get_cost(pre_d, s);
                    final_state = s;
                    break;
                }

                for (int i = 0; i < arraySize(_angles); i++) {
                    int theta_diff = _angles[i];
                    for (int v = 0; v < ops.size(); v++) { //check correctness of array size
                        if (_is_valid_op(theta_diff, ops[v])) {// TODO: translation is not always applicable
                            LocoState next_s = one_step.vr_op(ops[v], s, theta_diff, pre_d, _get_cost(pre_d, s));
                            int next_cost = one_step.get_cost();
                            int next_dist = one_step.get_dist();
                            if (_is_valid(next_s)) {
#ifdef DEBUG
                                std::cout << "----" << ops[v] << "\n";
                                s.print_state();
                                next_s.print_state(); printf("\n");
#endif // DEBUG
                                unsigned c = _get_cost(next_dist, next_s);
                                if (c > next_cost) {
                                    sigma[next_dist][next_s] = next_cost;
                                    parent_state[VRTuple(next_s, next_dist)] = VRTuple(s, pre_d);
                                    parent_op[VRTuple(next_s, next_dist)] = ops[v];
                                }
                            }
                        }
                    }
                }
            }

            if (complete_flag) {
                break;
            }
        }
        if (complete_flag) {
           _gen_vr_path(parent_state, parent_op, final_state);
        }
#ifdef DEBUG
        if (!complete_flag) {
            printf("[WARNNING] Not complete");
        }
#endif
    }
}

bool DPCvsp::_is_valid(LocoState s) {
    int gamma_p_x = s.getGammaPX();
    int gamma_p_y = s.getGammaPY();
    int gamma_v_x = s.getGammaVX();
    int gamma_v_y = s.getGammaVY();
    bool in_phy_space = (0 <= gamma_p_x &&  gamma_p_x < _phy_length && 0 <= gamma_p_y &&  gamma_p_y < _phy_width);
    bool in_vir_space = (0 <= gamma_v_x &&  gamma_v_x < _vir_length && 0 <= gamma_v_y &&  gamma_v_y < _vir_width);
    bool is_phy_obs = (_phy_point_map[Point(gamma_p_x, gamma_p_y)] == 0) ? true : false;
    bool is_vir_obs = (_vir_point_map[Point(gamma_v_x, gamma_v_y)] == 0) ? true : false;
    return in_phy_space && (!is_phy_obs) && in_vir_space && !is_vir_obs;
}

void DPCvsp::_gen_vr_path(std::map<VRTuple, VRTuple>& parent_state, std::map<VRTuple, std::string>& parent_op, LocoState s_t) {
    LocoState current_state = s_t;
    int dist = opt_path.distance;
    std::string pre_op = "null";
    while (dist >= 0) {
        opt_path.vr_path.push_back(current_state);
        opt_path.vr_ops.push_back(pre_op);
#ifdef DEBUG
        VRTuple t = parent_state[VRTuple(current_state, dist)];
        std::cout << "distance: " << dist << "\n";
        current_state.print_state();
        t.s.print_state();
#endif // DEBUG

        VRTuple pre_state = parent_state[VRTuple(current_state, dist)];
        pre_op = parent_op[VRTuple(current_state, dist)];
        dist = pre_state.d;
        current_state = pre_state.s;
    }
    std::reverse(opt_path.vr_path.begin(), opt_path.vr_path.end());
    std::reverse(opt_path.vr_ops.begin(), opt_path.vr_ops.end());
}

bool DPCvsp::_is_goal(LocoState s) {
    if (s.getGammaVX() == _gamma_t_v_x && s.getGammaVY() == _gamma_t_v_y) {
        return true;
    }
    return false;
}

bool DPCvsp::_is_valid_op(int theta_diff, std::string op) {
    bool valid = true;

    if (op == "translation") {
        if (theta_diff != 0 ) {
            valid = false;
        }
    }
    else if (op == "rotation") {
        if (theta_diff == 0) {
            valid = false;
        }
    }
    else if (op == "curvature") {
        if (theta_diff != 0) {
            valid = false;
        }
    }

    return valid;
}

unsigned DPCvsp::_get_cost(unsigned dist, LocoState s)
{
    auto it_1 = sigma.find(dist);
    if (it_1 == sigma.end())
    {
        return std::numeric_limits<unsigned>::max();
    }
    auto it_2 = it_1->second.find(s);
    if (it_2 == it_1->second.end())
    {
        return std::numeric_limits<unsigned>::max();
    }
    return it_2->second;
}


std::pair<float, float> DPCvsp::print_path() {
    std::vector<LocoState> path = opt_path.vr_path;
    std::vector<std::string> ops = opt_path.vr_ops;
    std::cout << opt_path.cost << "|";
    std::cout << opt_path.distance;
    for (int i = 0; i < path.size(); i++) {
        std::cout << "|" << path[i].to_string();
    }

    return std::make_pair(opt_path.distance, opt_path.cost);

}


DPCvsp::~DPCvsp()
{
}
