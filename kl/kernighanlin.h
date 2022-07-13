#ifndef COLAPREPROCESSOR_H
#define COLAPREPROCESSOR_H

#include <vector>
#include <tuple>
#include <map>
#include <set>

struct Point
{
    Point() : x(), y() {}
    Point(int a, int b) : x(a), y(b) {}
    int x;
    int y;
    friend bool operator<(const Point& l, const Point& r)
    {
        return std::tie(l.x, l.y) < std::tie(r.x, r.y);
    }
};

struct Edge
{
    Edge(unsigned p, unsigned l) : point(p), length(l) {}
    unsigned point;
    unsigned length;
    friend bool operator<(const Edge& l, const Edge& r)
    {
        return l.point < r.point;
    }
};

class KernighanLin
{
public:
    KernighanLin();

    void read_graph(const char* graph_file, size_t limit, float ratio);

private:

    void _read_graph(const char* graph_file);
    const Point& _get_point(unsigned idx) const;
    unsigned _get_idx(int x, int y) const;
    unsigned _get_idx(const Point& p) const;
    unsigned _get_idx_or_create_point(int x, int y);
    unsigned _get_idx_or_create_point(const Point& p);
    void _add_edge(unsigned s, unsigned t, unsigned length);

    void _write_gr(const char* file) const;
    void _write_co(const char* file) const;
    void _write_cola(const char* file, float ratio) const;
    void _write_partition(const char* file) const;
    void _write_cut_edge(const char* file) const;

    int _D(const std::set<unsigned>& A, const std::set<unsigned>& B, unsigned u) const;
    void _update_D(const std::set<unsigned>& A, const std::set<unsigned>& B, std::vector<int>& D) const;
    void _find_k_which_maximize_g_max(const std::vector<int>& gv, size_t &k_max, int &g_max) const;
    unsigned _c(unsigned a, unsigned b) const;
    void _kl(std::set<unsigned> &A0, std::set<unsigned> &B0) const;
    void _kl_recursive(std::set<unsigned> A, std::set<unsigned> B, size_t limit);
    void _kl_main(size_t limit);

    void _add_group(const std::set<unsigned> group);

    std::map<Point, unsigned> _point_map;
    std::vector<Point>        _points;
    std::vector< std::vector<Edge> > _graph;
    std::vector<unsigned>            _group_id;
    std::vector<std::set<unsigned>>  _groups;
};

#endif // COLAPREPROCESSOR_H
