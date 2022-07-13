#include "PhyIndex.h"
#include <math.h>
#include <limits>
#include <cmath>
#include <iomanip>
#include <list>

using namespace std;

#define _USE_MATH_DEFINES
#define SHAOHENG
//#define DEBUGGG

PhyIndex::PhyIndex()
{
    max_lv = 0;
    _angle_table = { {Step(1, 0), 0}, { Step(1, 1), 45}, {Step(0, 1), 90 },{ Step(-1, 1), 135 },
    { Step(-1, 0), 180 },{ Step(-1, -1), 225 },{ Step(0, -1), 270 },{ Step(1, -1), 315 } };

    phy_index_flag = false; // used by SH
}

PhyIndex::PhyIndex(int len)
{
    max_lv = len;
    _angle_table = { { Step(1, 0), 0 },{ Step(1, 1), 45 },{ Step(0, 1), 90 },{ Step(-1, 1), 135 },
    { Step(-1, 0), 180 },{ Step(-1, -1), 225 },{ Step(0, -1), 270 },{ Step(1, -1), 315 } };

    phy_index_flag = false; // used by SH
}

void PhyIndex::read_graph(const char* graph_file)
{
    string s(graph_file);
    if (s.empty())
    {
        return;
    }
    _read_graph(graph_file);
}

void PhyIndex::_read_graph(const char* graph_file)
{
    ifstream infile(graph_file);
    string line;
    unsigned type = 3;
    while (getline(infile, line))
    {
        istringstream iss(line);
        int x, y;
        if (!(iss >> x >> y)) {
            if (!(iss.str() == "pois" || iss.str() == "obs"
                || iss.str() == "width" || iss.str() == "length" || iss.str() == "refs")) {
                cout << "_read_graph::tags error" << "\n";
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
            getline(infile, line);
            istringstream iss(line);
            if (!(iss >> length)) {
                cout << "_read_graph::length error" << "\n";
            }
            continue;
        }
        else if (iss.str() == "width") {
            getline(infile, line);
            istringstream iss(line);
            if (!(iss >> width)) {
                cout << "_read_graph::width error" << "\n";
            }
            continue;
        }

        _insert_point(Point(x, y), type);
    }

    for (unsigned x = 0; x < length; x++) {
        for (unsigned y = 0; y < width; y++) {
            Point gamma_p = Point(x, y);
            if (_get_type(gamma_p) == 1000) { //is not obs, refs, and pois
                _point_map[gamma_p] = 3; //3: free space
            }
        }
    }
}

void PhyIndex::print_graph()
{
    for (map<Point, int>::iterator it = _point_map.begin(); it != _point_map.end(); ++it) {
        printf("(%d,%d), %d\n", it->first.x, it->first.y, it->second);
    }
}

int PhyIndex::_get_type(const Point& p) const
{
    auto ite = _point_map.find(p);
    if (ite == _point_map.end())
    {
        return  1000;
    }
    return ite->second;
}


void PhyIndex::_insert_point(const Point& p, unsigned type)
{
    _point_map.insert(make_pair(p, type));
}


void PhyIndex::export_sh_index(string out_file_name) {

    ofstream out_file;
    out_file.open(out_file_name);

    out_file << max_lv << "\n";
    out_file << "ab\n";
    for (int i = 1; i < max_lv + 1; ++i){
        out_file << _sh_alpha_labels[i] << " " << _sh_beta_labels[i] << "\n";
    }

    out_file << "Index\n";

    
    for (auto it = _sh_labels.begin(); it != _sh_labels.end(); ++it ){
        Point p1 = it->first.first;
        Point p2 = it->first.second;
        out_file << p1.x << " " << p1.y << " " << 0 << " " << p2.x << " " << p2.y << "\n";
        for (int i = 0; i < max_lv + 1; ++i){
            out_file << it->second[i] << " ";
        }
        out_file << "\n";        
    }

    out_file.close();

}

void PhyIndex::build_sh_index(bool simplify) {
    int infty = numeric_limits<int>::max();
    int table_size = max_lv + 1;

    //initialize length labels

    _sh_alpha_labels.push_back(0);
    _sh_beta_labels.push_back(0);

    for (int l=1; l<=max_lv; l++){
        _sh_alpha_labels.push_back( numeric_limits<float>::max() );
        _sh_beta_labels.push_back(0);
    }

    //initialize physical indexs
    for (auto it1 = _point_map.begin(); it1 != _point_map.end(); ++it1){ //starting point

        if (it1->second == 0) { continue; }
        Point s = it1->first;

        for (auto it2 = _point_map.begin(); it2 != _point_map.end(); ++it2){ //ending point

            if (it2->second == 0) { continue; }
            Point t = it2->first;

            vector<float> st_values;
            st_values.resize(table_size);
            _sh_labels[make_pair(s,t)] = st_values;            
        }
    }

    for(int i = 0; i < table_size; i++){
        for (auto it1 = _point_map.begin(); it1 != _point_map.end(); ++it1){ //starting point

            if (it1->second == 0) { continue; }
            Point s = it1->first;

            float tempbeta = infty;

            for (auto it2 = _point_map.begin(); it2 != _point_map.end(); ++it2){ //ending point

                if (it2->second == 0) { continue; }
                Point t = it2->first;

                float cost = Optimize_Step(s, t, 0, 0, i, 0).cost; // all thetas are set to 0; let's see how it works

                if (cost < _sh_alpha_labels[i]){
                    _sh_alpha_labels[i] = cost; // update alpha values
                }
                if (cost < tempbeta){
                    tempbeta = cost; // update tempbeta, which is min_t c(s,t,i)
                }
                _sh_labels[make_pair(s,t)][i] = cost; //store the label
            }

            if (tempbeta > _sh_beta_labels[i]){
                _sh_beta_labels[i] = tempbeta; // update beta values
            }
        }
    }

    //replace all infty values to -1
    for (auto it = _sh_labels.begin(); it != _sh_labels.end(); ++it){
        for (int i = 0; i < table_size; ++i){
            if (it->second[i] > 1000){
                it->second[i] = -1;
            }
        }
    }
}

map<Point, int>* PhyIndex::get_point_map()
{
    return &_point_map;
}

void PhyIndex::set_max_lv(int i){
    max_lv = i;
}

bool PhyIndex::is_valid_input(int x, int y)
{
    Point p(x,y);
    auto it = _point_map.find(p);
    if ( it == _point_map.end() ) { return false; }
    else {
        if ( it->second != 0 ) { return true; }
        else { return false; }
    }
}

float PhyIndex::RW_index_cost(Point s, Point t, float vir_steplength){

    if ((s.x == t.x) && (s.y == t.y)) {
        if (vir_steplength == 0) { return 0; }
        else { return numeric_limits<float>::max(); }
    } // if two states are equal

    float dx = t.x - s.x;
    float dy = t.y - s.y;
    float phy_steplength = pow( pow(dx,2)+pow(dy,2), 0.5);
    float trans_gain = float(vir_steplength)/phy_steplength;

    //cout << "trans_gain = " << trans_gain << endl;
    //cout << "temp Cost =" << RESET_COST + trans_cost(trans_gain) * vir_steplength << endl;
    return RESET_COST + trans_cost(trans_gain) * vir_steplength;
}

float PhyIndex::trans_cost(float gain){
    float hard_UB = TRANS_UB * 1.5;
    float hard_LB = TRANS_LB / 1.5;
    if ((gain > hard_UB) || (gain < hard_LB)){
        return 10000;
    }
    if (gain > TRANS_UB){
        return pow( (gain - TRANS_UB), 2);
    }
    else if (gain < TRANS_LB){
        return pow( (TRANS_LB - gain), 2);
    }
    else{
        return 0;
    }
}

float PhyIndex::rotate_cost(float gain){
    float hard_UB = ROTATE_UB * 1.5;
    float hard_LB = ROTATE_LB / 1.5;
    if ((gain > hard_UB) || (gain < hard_LB)){
        return 10000;
    }
    if (gain > ROTATE_UB){
        return pow( (gain - ROTATE_UB), 2);
    }
    else if (gain < ROTATE_LB){
        return pow( (ROTATE_LB - gain), 2);
    }
    else{
        return 0;
    }
}
/*
float PhyIndex::trans_cost(float gain){
    if (gain > TRANS_UB){
        return pow( (gain - TRANS_UB), 2);
    }
    else if (gain < TRANS_LB){
        return pow( (TRANS_LB - gain), 2);
    }
    else{
        return 0;
    }
}

float PhyIndex::rotate_cost(float gain){
    if (gain > ROTATE_UB){
        return pow( (gain - ROTATE_UB), 2);
    }
    else if (gain < ROTATE_LB){
        return pow( (ROTATE_LB - gain), 2);
    }
    else{
        return 0;
    }
}*/

pair<float,float> get_length_and_angle(Point p1, Point p2){
    float dx = p2.x - p1.x;
    float dy = p2.y - p1.y;
    float l = pow( pow(dx, 2.0) + pow(dy, 2.0), 0.5);
    float temp_angle = (atan2(dy, dx)/M_PI)*180.;
    float angle = (temp_angle<0)?temp_angle+360.:temp_angle;
    return make_pair(l, angle);
}

Point WalkOneStep(Point p, float l, float angle)
{
    angle *= M_PI;
    angle /= 180.;
    float x = p.x;
    float y = p.y;
    x += l * cos(angle);
    y += l * sin(angle);
    cout << "Phy: x = " << x << ", y = " << y << endl;
    p.x = round(x);
    p.y = round(y);
    return p;
}

FloatPoint WalkOneStep(FloatPoint p, float l, float angle)
{
    angle *= M_PI;
    angle /= 180.;
    p.x += l * cos(angle);
    p.y += l * sin(angle);
    cout << "Vir: x = " << p.x << ", y = " << p.y << endl;
    return p;
}

Point ToPoint(FloatPoint fp)
{
    Point p(round(fp.x), round(fp.y));
    return p;
}

float AddAngle(float a, float b)
{
    float s = a+b;
    if (s>360) { return s - 360.; }
    if (s<0) { return s + 360.; }
    return s;
}

RW PhyIndex::RW_cost(float vir_theta, float phy_theta, float vir_steplength, float vir_steptheta, float phy_steplength, float phy_steptheta, bool intermediate){
    if (phy_steplength == 0){
        RW rw = {100000.0, true, 1, 1};
        return rw;
    }

    float trans_gain = vir_steplength/phy_steplength;

    float cost_1 = RESET_COST + trans_cost(trans_gain) * vir_steplength;

    if(intermediate){
        RW rw = {cost_1, true, 1, trans_gain};
        return rw;
    }

    float vir_rotation = vir_steptheta - vir_theta;
    float phy_rotation = phy_steptheta - phy_theta;
    if (vir_rotation < 0) { vir_rotation += 360; }
    if (phy_rotation < 0) { phy_rotation += 360; }
    if (vir_rotation > 180 ) { // clockwise rotation
        vir_rotation = 360 - vir_rotation;
        phy_rotation = 360 - phy_rotation;
    }

    float cost_2 = 100000.0;
    float rotate_gain = 0;

    if (phy_rotation != 0){
        rotate_gain = vir_rotation/phy_rotation;
        cost_2 = rotate_cost(rotate_gain) + trans_cost(trans_gain) * vir_steplength;
    }

    if (cost_1 < cost_2){
        RW rw = {cost_1, true, 1, trans_gain};
        return rw;
    }
    else{
        RW rw = {cost_2, false, rotate_gain, trans_gain};
        return rw;
    }
}

Optimized_Step PhyIndex::Optimize_Step(Point s, Point t, float vir_theta, float phy_theta, float vir_steplength, float vir_steptheta){
    map<Point, vector<float> > costs;
    map<Point, vector<pair<Point, int> > > predecessors;

    int ceiling = ceil(vir_steplength);
    float scale = vir_steplength / float(ceiling);

    for (auto it = _point_map.begin(); it != _point_map.end(); ++it){
        if (it->second == 0) { continue; }
        else {
            Point p = it->first;
            vector<float> v1;
            v1.resize(ceiling + 1);
            v1[0] = 0;
            costs[p] = v1;
            vector<pair<Point, int> > v2;
            v2.resize(ceiling + 1);
            predecessors[p] = v2;
        }
    }

    for (auto it = _point_map.begin(); it != _point_map.end(); ++it){
        if (it->second == 0) { continue; }
        else {
            Point p = it->first;
            for(int i=1; i<=ceiling; i++){
                pair<float, float> la_s_p = get_length_and_angle(s, p);
                costs[p][i] = RW_cost(vir_theta, phy_theta, (scale * i), vir_steptheta, la_s_p.first, la_s_p.second, false).cost;
                predecessors[p][i] = make_pair(s, 0);
            }
        }
    }


    for (int i=1; i<=ceiling; i++){
        for (auto it = _point_map.begin(); it != _point_map.end(); ++it){
            if (it->second == 0) { continue; }
            else{
                Point p = it->first;
                for (auto it2 = _point_map.begin(); it2 != _point_map.end(); ++it2){
                    if (it2->second == 0) { continue; }
                    else{
                        Point q = it2->first;
                        for(int j=1; j<i; j++){
                            Point pred = predecessors[q][i-j].first;
                            pair<float, float> la_pred_q = get_length_and_angle(pred, q);
                            pair<float, float> la_p_q = get_length_and_angle(q, p);
                            float tempcost = costs[q][i-j] + RW_cost(vir_steptheta, la_pred_q.second, (scale * j),
                                vir_steptheta, la_p_q.first, la_p_q.second, true).cost;
                            if (tempcost < costs[p][i]){
                                costs[p][i] = tempcost;
                                predecessors[p][i] = make_pair(q, j);
                            }
                        }
                    }
                }
            }
        }
    }


    Point p = predecessors[t][ceiling].first;
    pair<float, float> la_p_t = get_length_and_angle(p, t);
    Optimized_Step opts = {costs[t][ceiling], la_p_t.second};
    return opts;
}

Optimized_Step PhyIndex::Optimize_Step(Point s, Point t, float vir_theta, float phy_theta, float vir_steplength, float vir_steptheta, 
                                       vector<pair<Point, vector<float> > > &temp_path){
    map<Point, vector<float> > costs;
    map<Point, vector<pair<Point, int> > > predecessors;
    map<Point, vector<vector<float> > > detail_op; // detail information of each step corresponding to costs

    //discretize
    int ceiling = ceil(vir_steplength);
    float scale = vir_steplength / float(ceiling);

    for (auto it = _point_map.begin(); it != _point_map.end(); ++it){ //init empty labels
        if (it->second == 0) { continue; }
        else {
            Point p = it->first;
            vector<float> v1;
            v1.resize(ceiling + 1);
            v1[0] = 0;
            costs[p] = v1;

            vector<vector<float> > v3;
            v3.resize(ceiling + 1);
            detail_op[p] = v3;
            for (auto it2 = detail_op[p].begin(); it2 != detail_op[p].end(); ++it2){
                it2->resize(3);
            }

            vector<pair<Point, int> > v2;
            v2.resize(ceiling + 1);
            predecessors[p] = v2;
        }
    }

    for (auto it = _point_map.begin(); it != _point_map.end(); ++it){ //init first step (allow rotation)
        if (it->second == 0) { continue; }
        else {
            Point p = it->first;
            for(int i=1; i<=ceiling; i++){
                pair<float, float> la_s_p = get_length_and_angle(s, p); //la_s_p is the (length, angle) pair of the vector s-p
                RW rw_p_i = RW_cost(vir_theta, phy_theta, (scale * i), vir_steptheta, la_s_p.first, la_s_p.second, false);
                costs[p][i] = rw_p_i.cost;
                // if reset => phy_steptheta; else -1;
                if (rw_p_i.reset == true) {
                    detail_op[p][i][0] = la_s_p.second;
                }
                else {
                    detail_op[p][i][0] = -1;
                }
                detail_op[p][i][1] = rw_p_i.rotate_gain;
                detail_op[p][i][2] = rw_p_i.translation_gain;
                predecessors[p][i] = make_pair(s, 0);
            }
        }
    }

    //DP
    for (int i=1; i<=ceiling; i++){
        for (auto it = _point_map.begin(); it != _point_map.end(); ++it){
            if (it->second == 0) { continue; }
            else{
                Point p = it->first;
                for (auto it2 = _point_map.begin(); it2 != _point_map.end(); ++it2){
                    if (it2->second == 0) { continue; }
                    else{
                        Point q = it2->first;
                        for(int j=1; j<i; j++){ //optimize the path to p via q
                            Point pred = predecessors[q][i-j].first;
                            pair<float, float> la_pred_q = get_length_and_angle(pred, q);
                            pair<float, float> la_p_q = get_length_and_angle(q, p);
                            RW temp_rw = RW_cost(vir_steptheta, la_pred_q.second, (scale * j), vir_steptheta, la_p_q.first, la_p_q.second, true);
                            float tempcost = costs[q][i-j] + temp_rw.cost;
                            
                            //since this is an intermediate operation, only the reset solution can be performed
                            if (tempcost < costs[p][i]){
                                costs[p][i] = tempcost;
                                //cout << "Cost to point (" << p.x << "," << p.y << ") in length " << i << " is updated to " << tempcost << "\n";
                                predecessors[p][i] = make_pair(q, i-j);
                                // if reset => phy_steptheta; else -1;
                                if (temp_rw.reset == true) {
                                    detail_op[p][i][0] = la_p_q.second;
                                }
                                else {
                                    detail_op[p][i][0] = -1;
                                }
                                detail_op[p][i][1] = temp_rw.rotate_gain;
                                detail_op[p][i][2] = temp_rw.translation_gain;
                            }
                        }
                    }
                }
            }
        }
    }

    //retrieving steps
    // TODO
    // predecessors <last point, the number of the step>
    vector<pair<Point, vector<float> > > phy_path;
    Point cur = t;
    phy_path.push_back( make_pair(t, detail_op[t][ceiling]) );
    Point last_point = predecessors[cur][ceiling].first;
    float last_shift = predecessors[cur][ceiling].second;
    do {
        //move then push
        cur = last_point;
        phy_path.push_back( make_pair(cur, detail_op[cur][last_shift]) );
        last_point = predecessors[cur][last_shift].first;
        last_shift = predecessors[cur][last_shift].second;
    } while (!(cur == s) || last_shift > 0);
    if ( !(cur == s) || last_shift != 0) { cout << "Wrong traversal, shift:" << last_shift << endl; }

    reverse(phy_path.begin(), phy_path.end());
    temp_path = phy_path;

    Point p = predecessors[t][ceiling].first;
    pair<float, float> la_p_t = get_length_and_angle(p, t);
    Optimized_Step opts = {costs[t][ceiling], la_p_t.second};
    return opts;
}

float PhyIndex::Greedy_Optimize_Path(Point s, float phy_theta, vector<pair<float, float> > vsteps, float vir_theta, Point vir_start, map<Point, float> vir_turning_angles){
    float totalcost = 0;
    Point now_point = s;
    float now_phy_theta = phy_theta;
    float now_vir_theta = vir_theta;

    vector<pair<float, float> > real_steps;

    for (auto it = vsteps.begin(); it != vsteps.end(); it++){
        int segments = ceil(it->first / max_lv);
        float seg_length = it->first / float(segments);
        //cout << "Step length = " << it->first << ", max_lv = " << max_lv << ", segments = " << segments << ", seglength = " << seg_length << endl;
        for (int i = 0; i < segments; i++){
            real_steps.push_back(make_pair(seg_length, it->second));
        }
    }

    // ^^^ partition vsteps into several real_steps

    vector<Point> ppath;
    ppath.push_back(s);
    vector<pair<Point, vector<float> > > phy_path; // the detailed operations of user study in physical space

    for (auto it = real_steps.begin(); it != real_steps.end(); it++){
        float stepcost = numeric_limits<float>::max();
        float phy_steptheta = 0;
        Point predecessor;
        float vir_steplength = it->first;
        float vir_steptheta = it->second;
        vector<pair<Point, vector<float> > > best_temp_path; //currecntly the best path
        for (auto it2 = _point_map.begin(); it2 != _point_map.end(); ++it2){
            if (it2->second == 0) { continue; }
            Point next_point = it2->first;
            vector<pair<Point, vector<float> > > temp_path;
            Optimized_Step opts = Optimize_Step(now_point, next_point, now_vir_theta, now_phy_theta, vir_steplength, vir_steptheta, temp_path);
            float tempcost = opts.cost;
            float temptheta = opts.theta;
            if (tempcost < stepcost){
                stepcost = tempcost;
                phy_steptheta = temptheta;
                predecessor = next_point;
                best_temp_path = temp_path;
            }
        }
        ppath.push_back(predecessor);
        totalcost += stepcost;
        now_phy_theta = phy_steptheta;
        now_vir_theta = vir_steptheta;
        now_point = predecessor;
        phy_path.insert( phy_path.end(), best_temp_path.begin(), best_temp_path.end() ); //update the best path to the result
    }


    //adjust and trim phy_path
    for (auto it =  phy_path.begin(); it != phy_path.end(); ++it){
        auto it2 = it;
        ++it2;
        Point coord1 = it->first;
        if ( it2 == phy_path.end() ){
            // clean the last content
            it->second[0] = -1.;
            it->second[1] = 1.;
            it->second[2] = 1.;
        }
        else {
            Point coord2 = it2->first;
            //cout << coord1.x << " " << coord1.y << " " << coord2.x << " " << coord2.y << endl;
            it->second[0] = it2->second[0];
            it->second[1] = it2->second[1];
            it->second[2] = it2->second[2];
        }
    }

    ofstream out_file;
    out_file.open(pathfilename, std::ofstream::app);
    cout << "Output to " << pathfilename << endl;

    //decide the center of physical play area, now hardcoded at (1,1)
    int xcenter = 1;
    int ycenter = 1;
    out_file << "Player Sart Pos: ";

    int xdiff = vir_start.x - s.x;
    int ydiff = vir_start.y - s.y;

    out_file << xcenter + xdiff << ", " << ycenter + ydiff << "\nPath:\n";

    //print the contents of phy_path
    cout << fixed << setprecision(4);
    for (auto it = phy_path.begin(); it != phy_path.end(); ++it){
        Point coord = it->first;
        vector<float> op = it->second;
        if ((op[0] == 0) && (op[1] == 0) && (op[2] == 0)){
            phy_path.erase(it);
            it--;
            continue;
        }
        cout << "(" << coord.x << "," << coord.y << ")\tRESET:" << op[0] << "\trotation:" << op[1] << "\ttranslation:" << op[2] << endl;
    }

    FloatPoint cur_v(vir_start.x, vir_start.y); //current virtual position
    Point cur_p = s; //current physical position
    float cur_phy_angle = 0.;
    float cur_vir_angle = 0.;
    float phy_angle_diff = 0.;
    float vir_angle_diff = 0.;
    float step_phy_angle = 0.;
    float step_vir_angle = 0.;
    float step_phy_length = 0.;
    float step_vir_length = 0.;

    for (auto it = phy_path.begin(); it != phy_path.end(); ++it){
        Point next_p;
        auto it2 = it;
        it2++;
        if (it2 != phy_path.end()){
            next_p = it2->first;
        }
        else{
            next_p = it->first;
        }
        vector<float> op = it->second;

        //if Reset used, but still rotation in the virtual world
        auto itv = vir_turning_angles.find(ToPoint(cur_v));
        float vir_ta = 0;
        if (itv != vir_turning_angles.end()) {
            vir_ta = itv->second; //will rotate even after reset
            cout << "vir_angle_diff (map) = " << vir_ta << endl;
        }

        pair<float,float> pla = get_length_and_angle(cur_p, next_p);
        step_phy_angle = pla.second;
        cout << "step_phy_angle = " << step_phy_angle << endl;
        step_phy_length = pla.first;

        if (op[0] != -1){ // reset
            op[0] = AddAngle(op[0], vir_ta * -1); // subtract the later rotation from reset
            vir_angle_diff = vir_ta;
            phy_angle_diff = AddAngle(op[0], cur_phy_angle * -1);
            phy_angle_diff = AddAngle(phy_angle_diff, vir_ta);
        }

        //step_phy_angle = AddAngle(step_phy_angle, vir_ta * -1); // subtract the later rotation from the physical target angle

        else{ // no reset
            phy_angle_diff = step_phy_angle - cur_phy_angle;
            if ((vir_ta != 0) && (vir_ta * phy_angle_diff < 0)){ // to avoid wrong direction of rotation
                if (phy_angle_diff < 0) {phy_angle_diff += 360.;}
                else if (phy_angle_diff > 0) { phy_angle_diff -= 360.;}
            }

            cout << "phy_angle_diff = " << phy_angle_diff << endl;
            vir_angle_diff = ((op[0] == -1) ? (phy_angle_diff * op[1]) : 0.);
            cout << "vir_angle_diff = " << vir_angle_diff << endl;
        }

        step_vir_length = step_phy_length * op[2];
        step_vir_angle = AddAngle(cur_vir_angle, vir_angle_diff);
        step_phy_angle = AddAngle(cur_phy_angle, phy_angle_diff);

        cout << "step_vir_angle = " << step_vir_angle << endl;
        cout << "step_phy_angle = " << step_phy_angle << endl;

        //for stop the rotation gain
        FloatPoint temp_v = WalkOneStep(cur_v, 0.5, step_vir_angle);

        cur_v = WalkOneStep(cur_v, step_vir_length, step_vir_angle);
        cur_p = WalkOneStep(cur_p, step_phy_length, step_phy_angle);
        cur_vir_angle = step_vir_angle;
        cur_phy_angle = step_phy_angle;

        if (!(cur_p == next_p)){
            cout << "Physical next point is wrong!!\n";
        }

        //output the RW gains
        out_file << "(X:" << cur_v.x << "/Y:1/Z:" << cur_v.y << "/Gt:" << op[2] << "/Gc:0/Gr:" << op[1] << "/Reset:" << op[0] << "),\n";

        //if a rotation gain is applied, stop the rotation gain once the user starts walking
        if (op[1] != 1){
            out_file << "(X:" << temp_v.x << "/Y:1/Z:" << temp_v.y << "/Gt:" << op[2] << "/Gc:0/Gr:" << 1 << "/Reset:" << -1 << "),\n";
        }

    }

    out_file.close();
    return totalcost;
}



PhyIndex::~PhyIndex()
{
}


