#include "MCPP.h"


MCPP::MCPP()
{
    _angle_table = { { Step(1, 0), 0 },{ Step(1, 1), 45 },{ Step(0, 1), 90 },{ Step(-1, 1), 135 },
    { Step(-1, 0), 180 },{ Step(-1, -1), 225 },{ Step(0, -1), 270 },{ Step(1, -1), 315 } };
}

std::vector<Step> MCPP::move_steps(VirPath vir_path){
    std::vector<Step> all_moves;
    for (int i = 1; i < vir_path.size();i++){
        std::vector<Step> a_move = _walk_straight(vir_path[i].x - vir_path[i-1].x, vir_path[i].y - vir_path[i-1].y);
        all_moves.insert(all_moves.end(), a_move.begin(), a_move.end());
    }

    return all_moves;
}

bool MCPP::is_valid_input(int x, int y)
{
    Point p(x,y);
    auto it = _point_map.find(p);
    if ( it == _point_map.end() ) { return false; }
    else {
        if ( it->second != 0 ) { return true; }
        else { return false; }
    }
}

void MCPP::read_graph(const char* graph_file)
{
    std::string s(graph_file);
    if (s.empty())
    {
	std::cout << "No graph file!\n";
        return;
    }
    _read_graph(graph_file);
}

void MCPP::_insert_point(const Point& p, unsigned type)
{
    _point_map.insert(std::make_pair(p, type));
}

int MCPP::_get_type(const Point& p) const
{
    auto ite = _point_map.find(p);
    if (ite == _point_map.end())
    {
        return  3;
    }
    return ite->second;
}

void MCPP::_read_graph(const char* graph_file)
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
            if (!(iss >> length)) {
                std::cout << "_read_graph::length error" << "\n";
            }
            continue;
        }
        else if (iss.str() == "width") {
            std::getline(infile, line);
            std::istringstream iss(line);
            if (!(iss >> width)) {
                std::cout << "_read_graph::width error" << "\n";
            }
            continue;
        }

        _insert_point(Point(x, y), type);
    }

    for (int x = 0; x < length; x++) {
        for (int y = 0; y < width; y++) {
            Point gamma_p = Point(x, y);
            if (_get_type(gamma_p) != 0) { //is not obs
                _insert_point(Point(x, y), 3); //3: free space
            }
        }
    }
}


unsigned MCPP::_get_angle(const Step& m) const
{
    auto ite = _angle_table.find(m);
    if (ite == _angle_table.end())
    {
        return std::numeric_limits<unsigned>::max();
    }
    return ite->second;
}

PhyPathes MCPP::find_min_cost_pathes(Point& p, std::vector<Step> vir_moves, int theta_p) {
    Pathes phy_min_cost_paths;
    Step one_move = vir_moves[0];
    Walk one_step = Walk();
    std::vector<std::string> ops = one_step.get_op();

    LocoState init_state = LocoState(0, 0, _get_angle(one_move), p.x, p.y, theta_p);
    std::map<Tuple, Tuple> parent_state;
    std::map<Tuple, std::string> parent_op;
    std::map<int, std::map<LocoState, unsigned>> sigma;

    sigma[0][init_state] = 0;
    parent_state[Tuple(init_state, 0)] = Tuple(init_state, -1);
    parent_op[Tuple(init_state, -1)] = "null";
    parent_state[Tuple(init_state, 0)].s.print_state(); printf("\n");
    bool complete_flag = false;
    for (int pre_d = 0; pre_d < vir_moves.size(); pre_d++) {
        std::map<LocoState, unsigned> current_states = sigma[pre_d];
        if (current_states.size()==0) {
            #ifdef DEBUG
            std::cout << "[WARNNING] no next states" << std::endl;
            #endif
            break;
        }
        for (auto it = current_states.begin(); it != current_states.end(); ++it)
        {
            LocoState s = it->first;

            int theta_diff;
            if (pre_d == 0) {
                theta_diff = 0;
            }
            else {
                theta_diff = (_get_angle(vir_moves[pre_d]) - _get_angle(vir_moves[pre_d-1]) + 360) % 360;
            }

            for (int v = 0; v < ops.size(); v++) {
                if (_is_valid_op(pre_d, vir_moves, ops[v])) {
                    LocoState next_s = one_step.vr_op(ops[v], s, theta_diff, pre_d, _get_cost(pre_d, s, sigma));
                    
                    int next_cost = one_step.get_cost();
                    int next_dist = one_step.get_dist();
                    if (_is_valid(next_s.getGammaPX(), next_s.getGammaPY())) {
    #ifdef DEBUG
                        std::cout << "----" << ops[v] << "\n";
                        s.print_state();
                        printf("\n Next:");
                        next_s.print_state();
                        printf("\n");
    #endif // DEBUG
                        unsigned c = _get_cost(next_dist, next_s, sigma);
                        
                        if ( c > next_cost) {
                            sigma[next_dist][next_s] = next_cost;                           
                            parent_state[Tuple(next_s, next_dist)] = Tuple(s, pre_d);
                            parent_op[Tuple(next_s, next_dist)] = ops[v];

                            if (next_dist == vir_moves.size()) {
                                complete_flag = true;
                            }
                        }
                        
                    }
                }
            }
            
        }
        
    }
    if (!complete_flag) {
        printf("[WARNNING] Not complete\n");
    }
    
    PhyPathes min_cost_path;
    if (complete_flag) {
        min_cost_path = _gen_phy_pathes(init_state, parent_state, parent_op, vir_moves.size(), sigma);
        return min_cost_path;
    }
    return min_cost_path;
}

unsigned MCPP::_get_cost(int dist, LocoState s, std::map<int, std::map<LocoState, unsigned>>& sigma)
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

bool MCPP::_is_valid(int gamma_p_x, int gamma_p_y) {
    bool in_space = (0 <= gamma_p_x &&  gamma_p_x < width && 0 <= gamma_p_y &&  gamma_p_y < length);
    bool is_obs = (_get_type(Point(gamma_p_x, gamma_p_y)) == 0) ? true : false;
    return in_space && !is_obs;
}

bool MCPP::_is_valid_op(int pre_d, std::vector<Step> vir_moves, std::string op) {
    bool valid = true;

    if (op == "translation 2") {//ensure two step straight paths
        int theta_diff_1 = -1;
        int theta_diff_2 = -1;
        if (pre_d == 0) {
            theta_diff_1 = 0;
        }
        else {
            theta_diff_1 = (_get_angle(vir_moves[pre_d]) - _get_angle(vir_moves[pre_d - 1]) + 360) % 360;
        }

        if (pre_d < vir_moves.size() - 1) {
            theta_diff_2 = (_get_angle(vir_moves[pre_d + 1]) - _get_angle(vir_moves[pre_d]) + 360) % 360;
        }

        if (theta_diff_1 != 0 || theta_diff_2 != 0 || pre_d == vir_moves.size() - 1) {
            valid = false;
        }
    }
    else if (op == "rotation 0.5") {
        int theta_diff;
        if (pre_d == 0) {
            theta_diff = 0;
        }
        else {
            theta_diff = (_get_angle(vir_moves[pre_d]) - _get_angle(vir_moves[pre_d - 1]) + 360) % 360;
        }

        if (theta_diff == 0 || (theta_diff % 90) != 0) {
            valid = false;
        }
    }
    else if (op == "curvature -45" || op == "curvature 45") {
        int theta_diff;
        if (pre_d == 0) {
            theta_diff = 0;
        }
        else {
            theta_diff = (_get_angle(vir_moves[pre_d]) - _get_angle(vir_moves[pre_d - 1]) + 360) % 360;
        }

        if (theta_diff != 0) {
            valid = false;
        }
    }

    return valid;
}

int MCPP::_get_min_cost(int dist, std::map<int, std::map<LocoState, unsigned>>& sigma) {
    int min_cost = std::numeric_limits<int>::max();;
    std::map<LocoState, unsigned> state_costs = sigma[dist];
    for (std::map<LocoState, unsigned>::iterator it = state_costs.begin(); it != state_costs.end(); ++it) {
        if (min_cost > it->second) {
            min_cost = it->second;
        }
    }
    return min_cost;
}

PhyPathes MCPP::_gen_phy_pathes(LocoState init_state, std::map<Tuple, Tuple> p_state, std::map<Tuple, std::string> p_op, 
    int vir_len, std::map<int, std::map<LocoState, unsigned>> sigma) {
    PhyPathes min_cost_pathes;
    if (sigma[vir_len].size() == 0) {
        return min_cost_pathes;
    }

    unsigned min_cost = _get_min_cost(vir_len, sigma);
    std::map<LocoState, unsigned> state_costs = sigma[vir_len];
    for (auto it = state_costs.begin(); it != state_costs.end(); ++it) {
        if (it->second == min_cost){
            PhyPath indexed_path;
            LocoState current_state = it->first;
            std::string pre_op = "null";
            int dist = vir_len;
            indexed_path.distance = dist;
            indexed_path.cost = 0;

            Tuple pre_state;
            while (dist >= 0) {
                indexed_path.physical_path.push_back(Point(current_state.getGammaPX(), current_state.getGammaPY()));
                indexed_path.vr_ops.push_back(pre_op);
                if (pre_op != "null") {
                    indexed_path.cost += RESET_180_COST;
                }
    #ifdef DEBUG
                Tuple t = parent_state[Tuple(current_state, dist)];
                std::cout << "distance: " << dist << "\n";
                current_state.print_state();
                t.s.print_state();
    #endif // DEBUG

                LocoState pre_s = p_state[Tuple(current_state, dist)].s;
                int d = p_state[Tuple(current_state, dist)].d;
                pre_op = p_op[Tuple(current_state, dist)];
                current_state = pre_s;
                dist = d;
            }

            std::reverse(indexed_path.physical_path.begin(), indexed_path.physical_path.end());
            std::reverse(indexed_path.vr_ops.begin(), indexed_path.vr_ops.end());

            min_cost_pathes.push_back(indexed_path);
        }
    }
    return min_cost_pathes;
}

std::vector<Step> MCPP::_walk_straight(int x, int y) {
    std::vector<Step> moves;
    int unit_vec_x = 0;
    int unit_vec_y = 0;

    if (x != 0 && y != 0) {
        unit_vec_x = int(x / std::abs(x));
        unit_vec_y = int(y / std::abs(y));
    }
    else if (x == 0 && y != 0) {
        unit_vec_x = 0;
        unit_vec_y = int(y / std::abs(y));
    }
    else if (y == 0 && x != 0) {
        unit_vec_x = int(x / std::abs(x));
        unit_vec_y = 0;
    }
    else {
        printf("PhyIndex ERROR: x_v - x_u == 0 and y_v - y_u == 0");
    }

    int x_i = 0;
    int y_i = 0;
    while (x_i != x && y_i != y) {
        moves.push_back(Step(unit_vec_x, unit_vec_y));
        x_i += unit_vec_x;
        y_i += unit_vec_y;
    }
    
    if (x_i == x && y_i != y){
        int delta_y = y - y_i;
        int d_y = (int)(delta_y / std::abs(delta_y));
        for (int y_j = y_i; std::abs(y_j) < std::abs(y); y_j += d_y) {
            moves.push_back(Step(0, d_y));
        }
    }
    else if (x_i != x && y_i == y) {
        int delta_x = x - x_i;
        int d_x = (int)(delta_x / std::abs(delta_x));
        for (int x_j = x_i; std::abs(x_j) < std::abs(x); x_j += d_x) {
            moves.push_back(Step(d_x, 0));
        }
    }

#ifdef DEBUG
    if (moves.size() != std::abs(x) && moves.size() != std::abs(y)) {
        printf("[ERROR] PhyIndex:_walk_straight wrong number of moves\n");
        printf("(%d, %d)\n", x, y);
        printf("size: %d\n", moves.size());
        exit(0);
    }
#endif
    
    return moves;
}

MCPP::~MCPP()
{
}
