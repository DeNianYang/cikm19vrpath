#include "VirIndex.h"

VirIndex::VirIndex(int m): max_step_length(m)
{
}

void VirIndex::read_graph(const char* graph_file)
{
    std::string s(graph_file);
    if (s.empty())
    {
        return;
    }
    _read_graph(graph_file);
}

void VirIndex::_read_graph(const char* graph_file)
{
    std::ifstream infile(graph_file);
    std::string line;
    unsigned type = 4;
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
            type = 3;
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
}

void VirIndex::print_graph()
{
    for (std::map<Point, int>::iterator it = _point_map.begin(); it != _point_map.end(); ++it) {
        printf("(%d,%d), %d\n", it->first.x, it->first.y, it->second);
    }
}

void VirIndex::_insert_point(const Point& p, unsigned type)
{
    _point_map.insert(std::make_pair(p, type));
}

void VirIndex::build_vir_index() {
    for (std::map<Point, int>::iterator it_1 = _point_map.begin(); it_1 != _point_map.end(); ++it_1)
    {
        for (std::map<Point, int>::iterator it_2 = _point_map.begin(); it_2 != _point_map.end(); ++it_2)
        {
            if ((it_1->second == 1 || it_1->second == 2) && (it_2->second == 1 || it_2->second == 2)) {// is poi or ref
                Point u = it_1->first;
                Point v = it_2->first;
                if ( (u.x != v.x || u.y != v.y)
                    && (visible_graph[u].find(v) == visible_graph[u].end()) ) {
                    if (_is_visible(u, v)) {
                        visible_graph[u][v] = Step(v.x - u.x, v.y - u.y);
                        visible_graph[v][u] = Step(u.x - v.x, u.y - v.y);
                    }
                }
            }
        }
    }
}

bool VirIndex::_is_visible(Point u, Point v) {

    int dx = v.x - u.x;
    int dy = v.y - u.y;
    float dist = std::pow(std::pow(dx,2)+std::pow(dy,2), 0.5);
    if (dist > max_step_length){
        return false;
    }

    bool visible = true;

    int unit_vec_x;
    int unit_vec_y;
    if ((v.x - u.x != 0) && (v.y - u.y != 0)) {
        unit_vec_x = int((v.x - u.x) / std::abs(v.x - u.x));
        unit_vec_y = int((v.y - u.y) / std::abs(v.y - u.y));
    }
    else if ((v.x - u.x == 0) && (v.y - u.y != 0)) {
        unit_vec_x = 0;
        unit_vec_y = int((v.y - u.y) / std::abs(v.y - u.y));
    }
    else if ((v.x - u.x != 0) && (v.y - u.y == 0)) {
        unit_vec_x = int((v.x - u.x) / std::abs(v.x - u.x));
        unit_vec_y = 0;
    }
    else {
        printf("[ERROR] VirIndex:_is_visible v.x - u.x == 0 and v.y - u.y == 0\n");
        exit(1);
    }


    int x_i = u.x;
    int y_i = u.y;
    std::vector<Point> diagonal_list;
    while (x_i != v.x && y_i != v.y) {
        diagonal_list.push_back(Point(x_i, y_i));
        x_i += unit_vec_x;
        y_i += unit_vec_y;
    }

    if (x_i == v.x && y_i != v.y) {
        int delta_y = v.y - y_i;
        int d_y = int((v.y - y_i) / std::abs(v.y - y_i));
        diagonal_list.push_back(Point(x_i, y_i));
        for (unsigned i = 0; i < diagonal_list.size(); i++) {
            int y_diff = 0;
            int x_d = diagonal_list[i].x;
            int y_d = diagonal_list[i].y;
            while (y_diff != delta_y + d_y) {
                if(_point_map[Point(x_d, y_d + y_diff)] == 3 ){
                    visible = false;
                    break;
                }
                y_diff += d_y;
            }
            if (!visible) {
                break;
            }
        }
    }
    else if (x_i != v.x && y_i == v.y) {
        int delta_x = v.x - x_i;
        int d_x = int((v.x - x_i) / std::abs(v.x - x_i));
        diagonal_list.push_back(Point(x_i, y_i));
        for (unsigned i = 0; i < diagonal_list.size(); i++) {
            int x_diff = 0;
            int x_d = diagonal_list[i].x;
            int y_d = diagonal_list[i].y;
            while (x_diff != (delta_x + d_x)) {
                if (_point_map[Point(x_d + x_diff, y_d)] == 3) {
                    visible = false;
                    break;
                }
                x_diff += d_x;
            }
            if (!visible) {
                break;
            }
        }
    }
    else if (x_i == v.x && y_i == v.y) {
        diagonal_list.push_back(Point(x_i, y_i));
        for (unsigned i = 0; i < diagonal_list.size(); i++) {
            int x_d = diagonal_list[i].x;
            int y_d = diagonal_list[i].y;
            if (_point_map[Point(x_d, y_d)] == 3) {
                visible = false;
                break;
            }
        }
    }

    return visible;
}

void VirIndex::export_to_file(std::string out_file_name) {
    std::ofstream out_file;
    out_file.open(out_file_name);

    for (std::map<Point, std::map<Point, Step>>::iterator it_1 = visible_graph.begin(); it_1 != visible_graph.end(); ++it_1) {
        Point u = it_1->first;
        for (std::map<Point, Step>::iterator it_2 = it_1->second.begin(); it_2 != it_1->second.end(); ++it_2) {
            Point v = it_2->first;
            Step vir_step = it_2->second;
            out_file << "---\n";
            out_file << u.x << " " << u.y << "\n";
            out_file << v.x << " " << v.y << "\n";
            out_file << vir_step.delta_x << " " << vir_step.delta_y << "\n";
        }
    }

    out_file.close();
}

void VirIndex::export_for_ksp(std::string out_file_name) {
    std::ofstream out_file;
    out_file.open(out_file_name);
    out_file << "{\n";

    for (auto it_1 = visible_graph.begin(); it_1 != visible_graph.end(); ++it_1) {
        Point u = it_1->first;
        if (it_1->second.begin() == it_1->second.end()) {
            continue;
        }
        out_file << "\""<< u.x << "," << u.y << "\"" << ": {";
        for (auto it_2 = it_1->second.begin(); it_2 != it_1->second.end(); ++it_2) {
            Point v = it_2->first;
            Step vir_step = it_2->second;
            float length = std::pow( std::pow(vir_step.delta_x, 2) + std::pow(vir_step.delta_y, 2), 0.5 );
            int l = std::ceil(length);
            if (std::next(it_2, 1) == it_1->second.end()) {
                out_file << "\"" << v.x << "," << v.y << "\"" << ": \"" << l << "\"}";
            }
            else {
                out_file << "\"" << v.x << "," << v.y << "\"" << ": \"" << l << "\",";
            }
        }

        if (std::next(it_1, 1) != visible_graph.end()) {
            out_file << ",\n";
        }
    }
    out_file << "\n}";
    out_file.close();
}

std::map<Point, int>* VirIndex::get_point_map()
{
    return &_point_map;
}

bool VirIndex::is_valid_input(int x, int y)
{
    Point p(x,y);
    auto it = _point_map.find(p);
    if ( it == _point_map.end() ) { return false; }
    else {
        if (it->second == 1) { return true; }
        else { return false; }
    }
}

bool VirIndex::is_obs(int x, int y)
{
    Point p(x,y);
    auto it = _point_map.find(p);
    if ( it == _point_map.end() ) { return false; }
    else {
        if (it->second == 3) { return true; }
        else { return false; }
    }
}

void VirIndex::output_user_study_map(std::string out_file_name)
{
    std::ofstream out_file;
    out_file.open(out_file_name);
    out_file << "Map Size: " << length << ',' << width << "\nMap:\n";
    for (int i = 0; i < length; i++){
        for (int j = 0; j < width; j++){
            out_file << (is_obs(i,j) ? "(1)" : "(0)");
            if (j != width-1){
                out_file << ", ";
            }
        }
        out_file << "\n";
    }
}

VirIndex::~VirIndex()
{
}
