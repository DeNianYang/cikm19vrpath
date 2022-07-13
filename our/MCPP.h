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

class MCPP
{
public:
    MCPP();
    ~MCPP();

    std::vector<Step> move_steps(VirPath vir_path);
    PhyPathes find_min_cost_pathes(Point& p, std::vector<Step> vir_moves, int theta_p);
    void read_graph(const char* graph_file);
    bool is_valid_input(int, int);

    unsigned width;
    unsigned length;


private:
    std::vector<Step> _walk_straight(int x, int y);
    PhyPathes _gen_phy_pathes(LocoState init_state, std::map<Tuple, Tuple> p_state, std::map<Tuple, std::string> p_op, int vir_len, std::map<int, std::map<LocoState, unsigned>> sigma);
    int _get_min_cost(int dist, std::map<int, std::map<LocoState, unsigned>>& sigma);
    unsigned _get_angle(const Step& m) const;
    bool _is_valid_op(int pre_d, std::vector<Step> vir_moves, std::string op);
    bool _is_valid(int gamma_p_x, int gamma_p_y);
    void _read_graph(const char* graph_file);
    void _insert_point(const Point& p, unsigned type);
    int _get_type(const Point& p) const;
    unsigned _get_cost(int dist, LocoState s, std::map<int, std::map<LocoState, unsigned>>& sigma);

    std::map<Step, int> _angle_table;
    std::map<Point, int> _point_map;
};

