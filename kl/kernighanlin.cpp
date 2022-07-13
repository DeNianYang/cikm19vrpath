#include "kernighanlin.h"
#include "include/rapidjson/filereadstream.h"
#include "include/rapidjson/filewritestream.h"
#include "include/rapidjson/document.h"
#include "include/rapidjson/error/en.h"
#include "include/rapidjson/writer.h"
#include <sstream>
#include <stdio.h>
#include <sys/time.h>
#include <algorithm>
#include <numeric>

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

double GetCurrentTimeSec()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

KernighanLin::KernighanLin() :
    _point_map(),
    _points(),
    _graph(),
    _group_id(),
    _groups()
{
    _get_idx_or_create_point(std::numeric_limits<int>::max(), std::numeric_limits<int>::max());
}

const Point& KernighanLin::_get_point(unsigned idx) const
{
    return _points.at(idx);
}

unsigned KernighanLin::_get_idx(int x, int y) const
{
    return _get_idx(Point(x, y));
}

unsigned KernighanLin::_get_idx(const Point& p) const
{
    auto ite = _point_map.find(p);
    if(ite == _point_map.end())
    {
        return std::numeric_limits<unsigned>::max();
    }
    return ite->second;
}

unsigned KernighanLin::_get_idx_or_create_point(int x, int y)
{
    return _get_idx_or_create_point(Point(x, y));
}

unsigned KernighanLin::_get_idx_or_create_point(const Point& p)
{
    auto ite = _point_map.find(p);
    if(ite != _point_map.end())
    {
        return ite->second;
    }

    unsigned id = _points.size();
    _point_map.insert(std::make_pair(p, id));
    _points.push_back(p);

    return id;
}

static Point _parse_point(const std::string& s)
{
    std::vector<std::string> strs = split(s, ',');
    int x = std::atoi(strs.at(0).c_str());
    int y = std::atoi(strs.at(1).c_str());

    return Point(x, y);
}

void KernighanLin::read_graph(const char* graph_file, size_t limit, float ratio)
{
    std::string s(graph_file);
    if(s.empty())
    {
        return;
    }

    _read_graph(graph_file);

    std::vector<std::string> strs = split(s, '.');

    _write_gr(  std::string(strs.at(0) + std::string(".gr")).c_str());
    _write_co(  std::string(strs.at(0) + std::string(".co")).c_str());
    _write_cola(std::string(strs.at(0) + std::string(".cola")).c_str(), ratio);

    _kl_main(limit);

    _write_partition(std::string(strs.at(0) + std::string(".partition")).c_str());
    _write_cut_edge( std::string(strs.at(0) + std::string(".cut_edge")).c_str());
}

static Point _read_point(const rapidjson::Value& v)
{
    assert(v.IsString());

    return _parse_point(v.GetString());
}

void KernighanLin::_read_graph(const char* graph_file)
{
    FILE* fp = fopen(graph_file, "r");

    if(nullptr == fp)
    {
        printf("Error: failed to open %s\n", graph_file);
        return;
    }

    char readBuffer[65536];
    rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
    rapidjson::Document d;
    if(d.ParseStream(is).HasParseError())
    {
        printf( "%s\n", rapidjson::GetParseError_En( d.GetParseError() ) );
        return;
    }

    _graph.clear();

    for ( rapidjson::Document::ConstMemberIterator ite = d.MemberBegin(); ite != d.MemberEnd(); ++ite )
    {
        // "4,1": {"6,0": "2", "4,2": "1", "3,0": "1"},
        const rapidjson::Value& node  = ite->name;
        const rapidjson::Value& edges = ite->value;
        Point p1   = _read_point(node);
        unsigned u = _get_idx_or_create_point(p1);
        (void) u;
        for ( rapidjson::Document::ConstMemberIterator edge_ite = edges.MemberBegin(); edge_ite != edges.MemberEnd(); ++edge_ite )
        {
            Point p2   = _read_point(edge_ite->name);
            unsigned v = _get_idx_or_create_point(p2);
            unsigned l = (unsigned) std::atoi(edge_ite->value.GetString());
            _add_edge(u, v, l);
        }
    }

    fclose(fp);

    for(std::vector<Edge>& edges : _graph)
    {
        std::sort(edges.begin(), edges.end());
    }
}

void KernighanLin::_add_edge(unsigned s, unsigned t, unsigned length)
{
    if(_graph.size() < s + 1)
    {
        _graph.resize(s + 1);
    }
    if(_graph.size() < t + 1)
    {
        _graph.resize(t + 1);
    }
    _graph[s].push_back(Edge(t, length));
}

// .gr file
// p sp n m (n: # of nodes, m: # of arcs)
// a U V W  (U, V: node id, W: arc weight )

void KernighanLin::_write_gr(const char* file) const
{
    FILE* fp = fopen(file, "w");

    if(nullptr == fp)
    {
        printf("Error: failed to open %s\n", file);
        return;
    }
    size_t n = _points.size() - 1;
    size_t m = 0;
    for(size_t i = 0; i < _graph.size(); ++i)
    {
        m += _graph[i].size();
    }
    fprintf(fp, "p sp %zu %zu\n", n, m);
    for(size_t i = 0; i < _graph.size(); ++i)
    {
        for(size_t j = 0; j < _graph[i].size(); ++j)
        {
            const Edge& edge = _graph[i][j];
            fprintf(fp, "a %zu %u %u\n", i, edge.point, edge.length);
        }
    }

    fclose(fp);
}

// .co file
// p aux sp co n (n: # of coordinate lines)
// v id X Y      (id: id of node, X,Y: X/Y-coordinates)

void KernighanLin::_write_co(const char* file) const
{
    FILE* fp = fopen(file, "w");

    if(nullptr == fp)
    {
        printf("Error: failed to open %s\n", file);
        return;
    }

    fprintf(fp, "p aux sp co %zu\n", _points.size() - 1);
    for(size_t i = 1; i < _points.size(); ++i)
    {
        fprintf(fp, "v %zu %d %d\n", i, _points[i].x, _points[i].y);
    }

    fclose(fp);
}

// node index start from zero
// # nodes
// u v1 weight cost  v2 weight cost ... -1
void KernighanLin::_write_cola(const char* file, float ratio) const
{
    FILE* fp = fopen(file, "w");

    if(nullptr == fp)
    {
        printf("Error: failed to open %s\n", file);
        return;
    }

    fprintf(fp, "%zu\n", _points.size() - 1);
    for(size_t i = 1; i < _graph.size(); ++i)
    {
        fprintf(fp, "%zu ", i);
        for(size_t j = 0; j < _graph[i].size(); ++j)
        {
            const Edge& edge = _graph[i][j];
	    float cost = float(edge.length) * ratio;
            int icost = std::ceil(cost);
            fprintf(fp, "%u %u %u ", edge.point, edge.length, icost);
        }
        fprintf(fp, "-1\n");
    }

    fclose(fp);
}

// there are total 2 partitions
// they are:
// ============================
// partition 0:
// 1 0.000000, 0.000001
// ...
// ============================
// partition 1:
// 11 0.000000, 0.000001
// ...
// ============================
void KernighanLin::_write_partition(const char* file) const
{
    FILE* fp = fopen(file, "w");

    if(nullptr == fp)
    {
        printf("Error: failed to open %s\n", file);
        return;
    }

    fprintf(fp, "there are total %zu partitions\n", _groups.size() - 1);
    fprintf(fp, "they are:\n");
    for(size_t i = 1; i < _groups.size(); ++i)
    {
        fprintf(fp, "============================\n");
        fprintf(fp, "partition %zu:\n", i - 1);
        const std::set<unsigned>& group = _groups[i];
        for(unsigned id : group)
        {
            fprintf(fp, "%u %u %u\n", id, _points[id].x, _points[id].y);
        }
    }
    fprintf(fp, "============================\n");

    fclose(fp);
}

void KernighanLin::_write_cut_edge(const char* file) const
{
    FILE* fp = fopen(file, "w");

    if(nullptr == fp)
    {
        printf("Error: failed to open %s\n", file);
        return;
    }

    for(size_t i = 1; i < _graph.size(); ++i)
    {
        const std::vector<Edge>& edges = _graph[i];
        for(const Edge& e : edges)
        {
            if(_group_id[i] != _group_id[e.point])
            {
                fprintf(fp, "%u %u\n", i, e.point);
            }
        }
    }

    fclose(fp);
}

int KernighanLin::_D(const std::set<unsigned>& A, const std::set<unsigned>& B, unsigned u) const
{
    assert(A.find(u) != A.end());
    int E = 0;
    int I = 0;
    for(const Edge& e : _graph.at(u))
    {
        assert(e.point != std::numeric_limits<unsigned>::max());
        if(A.find(e.point) != A.end())
        {
            I += e.length;
        }
        else if(B.find(e.point) != B.end())
        {
            E += e.length;
        }
    }
    return E - I;
}

void KernighanLin::_update_D(const std::set<unsigned>& A, const std::set<unsigned>& B, std::vector<int>& D) const
{
    std::fill(D.begin(), D.end(), 0);

    for(unsigned a : A)
    {
        D.at(a) = _D(A, B, a);
    }
    for(unsigned b : B)
    {
        D.at(b) = _D(B, A, b);
    }
}

void KernighanLin::_find_k_which_maximize_g_max(const std::vector<int>& gv, size_t& k_max, int& g_max) const
{
    k_max = std::numeric_limits<size_t>::max();
    g_max = std::numeric_limits<int>::min();
    int g = 0;
    for(size_t i = 0; i < gv.size(); ++i)
    {
        g += gv[i];
        if(g > g_max)
        {
            g_max = g;
            k_max = i;
        }
    }
}

unsigned KernighanLin::_c(unsigned a, unsigned b) const
{
    const std::vector<Edge>& edges = _graph.at(a);
    auto ite = std::lower_bound(edges.begin(), edges.end(), Edge(b, 0));
    if(ite != edges.end() && ite->point == b)
    {
        return ite->length;
    }
    return 0;
}

void KernighanLin::_kl(std::set<unsigned>& A0, std::set<unsigned>& B0) const
{
    while(true)
    {
        std::set<unsigned> A = A0;
        std::set<unsigned> B = B0;
        std::vector<int> D(_points.size(), 0);
        _update_D(A, B, D);

        std::vector<int> gv;
        std::vector<unsigned> av;
        std::vector<unsigned> bv;

        size_t size = A.size() > B.size() ? B.size() : A.size();
        for(size_t i = 0; i < size; ++i)
        {
            int max_g      = std::numeric_limits<int>::min();
            unsigned max_a = std::numeric_limits<unsigned>::max();
            unsigned max_b = std::numeric_limits<unsigned>::max();
            for(unsigned a : A)
            {
                for(unsigned b : B)
                {
                    int g = D.at(a) + D.at(b) - 2 * _c(a, b);
                    if(g > max_g)
                    {
                        max_g = g;
                        max_a = a;
                        max_b = b;
                    }
                }
            }
            if(std::numeric_limits<int>::min() == max_g)
            {
                continue;
            }
            gv.push_back(max_g);
            av.push_back(max_a);
            bv.push_back(max_b);
            A.erase(max_a);
            B.erase(max_b);
            _update_D(A, B, D);
        }

        size_t k_max = std::numeric_limits<size_t>::max();
        int g_max    = std::numeric_limits<int>::min();
        _find_k_which_maximize_g_max(gv, k_max, g_max);

        if(g_max <= 0)
        {
            break;
        }

        for(size_t i = 0; i <= k_max; ++i)
        {
            A0.erase(av.at(i));
            B0.erase(bv.at(i));
            A0.insert(bv.at(i));
            B0.insert(av.at(i));
        }
    }
}

void KernighanLin::_kl_recursive(std::set<unsigned> A, std::set<unsigned> B, size_t limit)
{
    _kl( A, B);

    if(A.size() <= limit && B.size() <= limit)
    {
        _add_group(A);
        _add_group(B);
        return;
    }

    {
        std::vector<unsigned> v(A.begin(), A.end());
        std::set<unsigned> A0(v.begin(), v.begin() + v.size() / 2);
        std::set<unsigned> B0(v.begin() + v.size() / 2, v.end());
        _kl_recursive( A0, B0, limit);
    }
    {
        std::vector<unsigned> v(B.begin(), B.end());
        std::set<unsigned> A0(v.begin(), v.begin() + v.size() / 2);
        std::set<unsigned> B0(v.begin() + v.size() / 2, v.end());
        _kl_recursive( A0, B0, limit);
    }
}

void KernighanLin::_kl_main(size_t limit)
{
    _group_id.clear();
    _group_id.resize(_points.size(), std::numeric_limits<unsigned>::max());

    _add_group(std::set<unsigned>({0u}));

    std::vector<unsigned> index(_points.size() - 1, 0);
    std::iota(index.begin(), index.end(), 1);
    std::set<unsigned> A(index.begin(), index.begin() + index.size() / 2);
    std::set<unsigned> B(index.begin() + index.size() / 2, index.end());

    _kl_recursive( A, B, limit);
}

void KernighanLin::_add_group(const std::set<unsigned> group)
{
    unsigned id = _groups.size();
    for(unsigned a : group)
    {
        _group_id[a] = id;
    }
    _groups.push_back(group);
}
