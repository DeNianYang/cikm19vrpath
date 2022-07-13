#ifndef PhyIndex_H
#define PhyIndex_H

#include <tuple>

#include "Tools.h"
#include "Walk.h"
#include "LocoState.h"


struct Tuple
{
    Tuple() : s(), d() {}
    Tuple(LocoState a, int b) : s(a), d(b) {}
    LocoState s;
    int d;
    friend bool operator==(const Tuple& l, const Tuple& r)
    {
        return std::tie(l.s, l.d) == std::tie(r.s, r.d);
    }

    friend bool operator<(const Tuple& l, const Tuple& r)
    {
        return std::tie(l.s, l.d) < std::tie(r.s, r.d);
    }
};

class PhyIndex
{
public:
    PhyIndex();
    PhyIndex(int vir_len);
    ~PhyIndex();
    void read_graph(const char* graph_file);
    void build_phy_index();
    void print_graph();
    void export_to_file(std::string out_file);

    //Shao-Heng's index functions
    void build_sh_index(bool);
    //void trim_sh_index();
    void export_sh_index(std::string out_file);
    void set_max_lv(int);
    float RW_index_cost(Point, Point, float);
    float trans_cost(float);
    float rotate_cost(float);
    RW RW_cost(float vir_theta, float phy_theta,
    float vir_steplength, float vir_steptheta, float phy_steplength, float phy_steptheta, bool intermediate = false);
    Optimized_Step Optimize_Step(Point s, Point t, float vir_theta, float phy_theta, float vir_steplength, float vir_steptheta);
    Optimized_Step Optimize_Step(Point s, Point t, float vir_theta, float phy_theta, float vir_steplength, float vir_steptheta, std::vector<std::pair<Point, 
                                 std::vector<float> > > &temp_path);
    float Greedy_Optimize_Path(Point s, float phy_theta, std::vector<std::pair<float, float> > vsteps, float vir_theta, Point vir_start, std::map<Point, float> vir_turning_angles);

    bool is_valid_input(int, int);

    std::map<Point, int>* get_point_map();

    friend class SearchSpace;
    unsigned width;
    unsigned length;
    std::string pathfilename;

private:
    void _read_graph(const char* graph_file);
    int _get_type(const Point& p) const;
    void _insert_point(const Point& p, unsigned type);
    /*
    void _insert_label_set(const Point& p, const PhyLabel& l);
    PhyLabel _build_physical_label_set(Point& p, int theta_p);
    std::vector<Step> _walk_straight(int x, int y);
    PhyPathes _MCPP(Point& p, std::vector<Step> vir_moves, int theta_p);
    unsigned _get_angle(const Step& m) const;
    unsigned _get_cost(int dist, LocoState s, std::map<int, std::map<LocoState, unsigned>>& sigma);
    bool _is_valid(int gamma_p_x, int gamma_p_y);
    PhyPathes _gen_phy_pathes(LocoState init_state, std::map<Tuple, Tuple> parent_state, std::map<Tuple, std::string> parent_op, int vir_len, std::map<int, std::map<LocoState, unsigned>> sigma);
    int _get_min_cost(int dist, std::map<int, std::map<LocoState, unsigned>>& sigma);
    bool _is_valid_op(unsigned pre_d, std::vector<Step> vir_moves, std::string op);
    */


    int max_lv;

    std::map<Point, int> _point_map;
    //std::map<PhyState, PhyLabel> _phy_labels;
    std::map<Step, int> _angle_table;

    //Shao-Heng's index
    std::map<std::pair<Point, Point> , std::vector<float> > _sh_labels;
    std::vector<float> _sh_alpha_labels;
    std::vector<float> _sh_beta_labels;


    const PhyState* getState(int x, int y, int angle);
    bool phy_index_flag;
};

#endif // PhyIndex_H
