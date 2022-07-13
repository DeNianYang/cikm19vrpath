#ifndef VirIndex_H
#define VirIndex_H

#include "Tools.h"

class VirIndex
{
public:
    //VirIndex();
    VirIndex(int);
    ~VirIndex();
    void read_graph(const char* graph_file);
    void build_vir_index();
    void export_to_file(std::string out_file_name);
    void export_for_ksp(std::string out_file_name);
    void print_graph();
    void output_user_study_map(std::string out_file_name);

    std::map<Point, int>* get_point_map();

    bool is_obs(int, int);
    bool is_valid_input(int, int);
    unsigned width;
    unsigned length;

private:
    void _read_graph(const char* graph_file);
    void _insert_point(const Point& p, unsigned type);
    bool _is_visible(Point u, Point v);


    std::map<Point, int> _point_map;
    std::map<Point, std::map<Point, Step>> visible_graph;

    int max_step_length;
};
#endif

