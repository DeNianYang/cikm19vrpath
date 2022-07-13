#ifndef KSHORTESTPATH_H
#define KSHORTESTPATH_H

#include "Tools.h"

class KShortestPath
{
public:
    KShortestPath();

    void read_graph(const char* graph_file, unsigned K);
    bool construct_index();
    VirPathes query(int x1, int y1, int x2, int y2) const;
    void print_summary() const;
    void read_index(const char* file);
    void write_index(const char* file) const;
    void add_indexing_time(double t);
    void add_query_time(double t);
    void read_type(const char* map_file);

    // test use
    void generate_sample_graph();
    void print_member_data() const;
    VirPathes query(unsigned s, unsigned t) const;

    double get_indexing_time() const { return _indexing_time; }
    double get_loop_count_time() const { return _loop_count_time; }
    double get_query_time() const { return _query_time; }
    size_t get_memory_usage() const;

private:

    void _read_graph(const char* graph_file);
    //    Point _parse_point(const std::string& s) const;
    const Point& _get_point(unsigned idx) const;
    unsigned _get_idx(int x, int y) const;
    unsigned _get_idx(const Point& p) const;
    unsigned _get_idx_or_create_point(int x, int y);
    unsigned _get_idx_or_create_point(const Point& p);
    bool _add_point_and_edges(const std::vector<int>& inputs);
    void _add_edge(unsigned s, unsigned t, unsigned length);
    void _add_double_edge(unsigned s, unsigned t, unsigned length);

    void _init();
    bool _labeling();
    void _count_loop(unsigned s, bool& status);
    void _pruned_bfs(unsigned s, bool& status);
    kPathes _get_pathes(const kPathes& pathes_s, const kPathes& pathes_v, const kPathes& pathes_t) const;
    DistancePathes _query(unsigned s, unsigned t) const;
    bool _fast_pruning(unsigned s, unsigned v, unsigned distance) const;
    bool _pruning(unsigned s, unsigned v, unsigned distance) const;
    void _read_type(const char* map_file);
    void _insert_map(const Point& p, unsigned type);

    VirPath _convert_path_point(const kPath& path) const;
    void _print_path(const kPath& path) const;
    void _print_label(const kLabel& label) const;

    size_t _get_path_memory(const kPath& path) const;
    size_t _get_pathes_memory(const kPathes& pathes) const;
    size_t _get_label_memory(const kLabel& label) const;


    std::map<Point, unsigned> _point_map;
    std::vector<Point>        _points;
    std::vector< std::vector<Edge> > _graph;
    std::vector< kLabel >                      _loop_label;   // v -> v : distance, pathes ...
    std::vector< std::map<unsigned, kLabel > > _dist_label;   // v -> u : distance, pathes ...

    size_t _K; // num of shortest path
    size_t _V; // num of vertex

    double _indexing_time;
    double _loop_count_time;
    mutable double _query_time;
    std::map<Point, int> _type_map;

    static const unsigned INF;
};

#endif // KSHORTESTPATH_H
