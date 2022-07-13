#include "kshortestpath.h"
#include "include/rapidjson/filereadstream.h"
#include "include/rapidjson/filewritestream.h"
#include "include/rapidjson/document.h"
#include "include/rapidjson/error/en.h"
#include "include/rapidjson/writer.h"
#include <fstream>
#include <sstream>
#include <limits>
#include <queue>
#include <stdio.h>
#include <algorithm>
#include <iterator>
#include <ctime>

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

const unsigned KShortestPath::INF = std::numeric_limits<unsigned>::max() / 2;
//const unsigned KShortestPath::INF = 128;

KShortestPath::KShortestPath() :
    _point_map(),
    _points(),
    _graph(),
    _loop_label(),
    _dist_label(),
    _K(),
    _V(),
    _indexing_time(0.0),
    _loop_count_time(0.0),
    _query_time(0.0)
{
}

const Point& KShortestPath::_get_point(unsigned idx) const
{
    return _points.at(idx);
}

unsigned KShortestPath::_get_idx(int x, int y) const
{
    return _get_idx(Point(x, y));
}

unsigned KShortestPath::_get_idx(const Point& p) const
{
    auto ite = _point_map.find(p);
    if (ite == _point_map.end())
    {
        return std::numeric_limits<unsigned>::max();
    }
    return ite->second;
}

unsigned KShortestPath::_get_idx_or_create_point(int x, int y)
{
    return _get_idx_or_create_point(Point(x, y));
}

unsigned KShortestPath::_get_idx_or_create_point(const Point& p)
{
    auto ite = _point_map.find(p);
    if (ite != _point_map.end())
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

void KShortestPath::read_graph(const char* graph_file, unsigned K)
{
    _K = K;

    _read_graph(graph_file);
}

static Point _read_point(const rapidjson::Value& v)
{
    assert(v.IsString());

    return _parse_point(v.GetString());
}

void KShortestPath::_read_graph(const char* graph_file)
{
    FILE* fp = fopen(graph_file, "r");

    char readBuffer[65536];
    rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
    rapidjson::Document d;
    if (d.ParseStream(is).HasParseError())
    {
        printf("%s\n", rapidjson::GetParseError_En(d.GetParseError()));
        return;
    }

    _graph.clear();

    for (rapidjson::Document::ConstMemberIterator ite = d.MemberBegin(); ite != d.MemberEnd(); ++ite)
    {
        // "4,1": {"6,0": "2", "4,2": "1", "3,0": "1"},
        const rapidjson::Value& node = ite->name;
        const rapidjson::Value& edges = ite->value;
        Point p1 = _read_point(node);
        unsigned u = _get_idx_or_create_point(p1);
        (void)u;
        for (rapidjson::Document::ConstMemberIterator edge_ite = edges.MemberBegin(); edge_ite != edges.MemberEnd(); ++edge_ite)
        {
            Point p2 = _read_point(edge_ite->name);
            unsigned v = _get_idx_or_create_point(p2);
            unsigned l = (unsigned)std::atoi(edge_ite->value.GetString());
            _add_edge(u, v, l);
        }
    }

    fclose(fp);
}

void KShortestPath::_add_edge(unsigned s, unsigned t, unsigned length)
{
    if (_graph.size() < s + 1)
    {
        _graph.resize(s + 1);
    }
    if (_graph.size() < t + 1)
    {
        _graph.resize(t + 1);
    }
    _graph[s].push_back(Edge(t, length));
}

void KShortestPath::_add_double_edge(unsigned s, unsigned t, unsigned length)
{
    _add_edge(s, t, length);
    _add_edge(t, s, length);
}

void KShortestPath::generate_sample_graph()
{
    _graph.clear();

    // example.txt

    _K = 16;

    _add_double_edge(1, 3, 1);
    _add_double_edge(1, 6, 1);
    _add_double_edge(2, 3, 1);
    _add_double_edge(2, 4, 1);
    _add_double_edge(2, 8, 1);
    _add_double_edge(2, 12, 1);
    _add_double_edge(3, 6, 1);
    _add_double_edge(3, 7, 1);
    _add_double_edge(3, 11, 1);
    _add_double_edge(4, 12, 1);
    _add_double_edge(5, 10, 1);
    _add_double_edge(6, 10, 1);
    _add_double_edge(8, 9, 1);
    _add_double_edge(8, 12, 1);
    _add_double_edge(9, 10, 1);
    _add_double_edge(11, 9, 1);
    _add_double_edge(12, 9, 1);
}

VirPath KShortestPath::_convert_path_point(const kPath& path) const
{
    VirPath point_path;
    for (size_t j = 0; j < path.size(); ++j)
    {
        const Point& p = _get_point(path[j]);
        point_path.push_back(p);
    }

    return point_path;
}

void KShortestPath::_print_path(const kPath& path) const
{
    printf("    ");
    for (size_t j = 0; j < path.size(); ++j)
    {
        printf("%u ", path[j]);
        if (j < path.size() - 1)
        {
            printf("-> ");
        }
        else
        {
            printf("\n");
        }
    }
}

void KShortestPath::_print_label(const kLabel& label) const
{
    for (const auto& kv : label)
    {
        printf("  (distance %u, count %zu):\n", kv.first, kv.second.size());
        const kPathes& pathes = kv.second;
        for (const kPath& path : pathes)
        {
            _print_path(path);
        }
    }
}

void KShortestPath::print_member_data() const
{
    for (size_t i = 0; i < _V; ++i)
    {
        printf("node %zu: \n", i);
        printf("  loop label:\n");
        _print_label(_loop_label[i]);

        printf("  dist label:\n");
        const std::map<unsigned, kLabel >& map = _dist_label[i];
        for (const auto& kv : map)
        {
            printf("    %zu -> %u\n", i, kv.first);
            _print_label(kv.second);
        }
    }
}

void KShortestPath::print_summary() const
{
    printf("Loop count time : %f s\n", get_loop_count_time());
    printf("Indexing time   : %f s\n", get_indexing_time());
    printf("Query time      : %f s\n", get_query_time());
    size_t mem_size = get_memory_usage();
    if (mem_size >> 20)
    {
        printf("Memory          : %zu MByte (%zu Byte)\n", mem_size >> 20, mem_size);
    }
    else if (mem_size >> 10)
    {
        printf("Memory          : %zu KByte (%zu Byte)\n", mem_size >> 10, mem_size);
    }
    else
    {
        printf("Memory          : %zu Byte\n", mem_size);
    }
}

VirPathes KShortestPath::query(unsigned s, unsigned t) const
{
    VirPathes point_pathes;
    if (s >= _V || t >= _V)
    {
        return point_pathes;
    }
    const Point& p1 = _get_point(s);
    const Point& p2 = _get_point(t);
    printf("query (%d, %d) -> (%d, %d)\n", p1.x, p1.y, p2.x, p2.y);

    _query_time = -GetCurrentTimeSec();
    DistancePathes pathes = _query(s, t);
    _query_time += GetCurrentTimeSec();

    for (const DistancePath& d_path : pathes)
    {
        printf("distance: %u path: ", d_path.distance);
        const kPath& path = d_path.path;
        VirPath point_path = _convert_path_point(path);
        point_pathes.push_back(point_path);
    }

    return point_pathes;
}

VirPathes KShortestPath::query(int x1, int y1, int x2, int y2) const
{
    unsigned s = _get_idx(x1, y1);
    unsigned t = _get_idx(x2, y2);
    VirPathes pathes = query(s, t);
    return pathes;
}

bool KShortestPath::construct_index()
{
    _init();
    return _labeling();
}

void KShortestPath::_init()
{
    _V = _graph.size();

    _loop_label.resize(_V);
    _dist_label.resize(_V);
    for (size_t i = 0; i < _dist_label.size(); ++i)
    {
        // add v ->v : distance 0
        _dist_label[i][i][0].push_back(kPath(1, i));
    }
}

bool KShortestPath::_labeling()
{
    bool status = true;
    _loop_count_time = -GetCurrentTimeSec();
    for (size_t i = 0; i < _V; ++i)
    {
        if (_type_map[_get_point(i)] == 1) {
            //std::cout << "Counting loop for point " << i << std::endl;
            _count_loop(i, status);
        }
        //_count_loop(i, status);
    }
    _loop_count_time += GetCurrentTimeSec();

    _indexing_time = -GetCurrentTimeSec();
    for (size_t i = 0; i < _V; ++i)
    {
        if (_type_map[_get_point(i)] == 1) {
            //std::cout << "Doing pruned BFS for point " << i << std::endl;
            _pruned_bfs(i, status);
        }
        //_pruned_bfs(i, status);
    }
    _indexing_time += GetCurrentTimeSec();

    return status;
}

void KShortestPath::_count_loop(unsigned s, bool& status)
{
    size_t count = 0;

    std::priority_queue<DistancePath> node_que;
    node_que.push(DistancePath(s));

    clock_t tt = clock();

    while (!node_que.empty() && count < _K)
    {
        clock_t nowt = clock();
        if ( float(nowt - tt)/CLOCKS_PER_SEC > 60 ){
            printf("Warning: Counting loop too long. Dangerous. Break.\n");
            status = false;
            break;
        }

        DistancePath dp = node_que.top(); node_que.pop();
        unsigned v = dp.path.back();
        if (dp.distance == INF)
        {
            if (status)
            {
                printf("Warning: Self loops become too long.\n");
            }
            status = false;
            break;
        }

        if (v == s)
        {
            _loop_label[s][dp.distance].push_back(dp.path);
            ++count;
        }

        for (size_t i = 0; i < _graph[v].size(); ++i)
        {
            const Edge& e = _graph[v][i];
            if (e.point < s)
            {
                continue;
            }
            DistancePath d(dp);
            d.distance += e.length;
            d.path.push_back(e.point);
            node_que.push(d);
        }
    }
}

void KShortestPath::_pruned_bfs(unsigned s, bool& status)
{
    size_t count = 0;

    std::priority_queue<DistancePath> node_que;
    node_que.push(DistancePath(s));

    std::vector<bool> pruned(_V, false);
    std::map<unsigned, kLabel >& map_s = _dist_label.at(s);

    while (!node_que.empty() && count < _K)
    {
        DistancePath dp = node_que.top(); node_que.pop();
        unsigned v = dp.path.back();

        if (pruned.at(v))
        {
            continue;
        }

        if (dp.distance == INF)
        {
            if (status)
            {
                printf("Warning: Distance from a source node becomes too long.\n");
            }
            status = false;
            break;
        }

        pruned.at(v) = _pruning(s, v, dp.distance);
        if (pruned.at(v))
        {
            continue;
        }

        map_s[v][dp.distance].push_back(dp.path);
        std::map<unsigned, kLabel >& map_v = _dist_label.at(v);
        map_v[s][dp.distance].push_back(kPath(dp.path.rbegin(), dp.path.rend()));

        for (size_t i = 0; i < _graph[v].size(); ++i)
        {
            const Edge& e = _graph[v][i];
            if (e.point < s)
            {
                continue;
            }
            DistancePath d(dp);
            d.distance += e.length;
            d.path.push_back(e.point);
            node_que.push(d);
        }
    }
}

kPathes KShortestPath::_get_pathes(const kPathes& pathes_s, const kPathes& pathes_v, const kPathes& pathes_t) const
{
    kPathes pathes;
    for (const kPath& path_s : pathes_s)
    {
        for (const kPath& path_v : pathes_v)
        {
            for (const kPath& path_t : pathes_t)
            {
                kPath path = path_s;
                if (path_v.size() > 1)
                {
                    path.insert(path.end(), path_v.begin() + 1, path_v.end());
                }
                if (path_t.size() > 1)
                {
                    if (path_t.front() == path.back())
                    {
                        path.insert(path.end(), path_t.begin() + 1, path_t.end());
                    }
                    else
                    {
                        path.insert(path.end(), path_t.rbegin() + 1, path_t.rend());
                    }
                }
                pathes.push_back(path);
            }
        }
    }

    return pathes;
}

DistancePathes KShortestPath::_query(unsigned s, unsigned t) const
{
    const std::map<unsigned, kLabel>& map_s = _dist_label.at(s);
    const std::map<unsigned, kLabel>& map_t = _dist_label.at(t);

    kLabel result;

    for (const auto& kv : map_s)
    {
        auto ite = map_t.find(kv.first);
        if (ite == map_t.end())
        {
            continue;
        }
        unsigned v = kv.first;
        const kLabel& label_s = kv.second;
        const kLabel& label_t = ite->second;
        const kLabel& label_v = _loop_label.at(v);

        for (const auto& kv_s : label_s)
        {
            for (const auto& kv_v : label_v)
            {
                for (const auto& kv_t : label_t)
                {
                    unsigned distance = kv_s.first + kv_v.first + kv_t.first;
                    kPathes  pathes = _get_pathes(kv_s.second, kv_v.second, kv_t.second);
                    kPathes& ref_pathes = result[distance];
                    ref_pathes.insert(ref_pathes.end(), pathes.begin(), pathes.end());
                }
            }
        }
    }

    DistancePathes pathes;
    for (const auto& kv : result)
    {
        for (const kPath& path : kv.second)
        {
            pathes.push_back(DistancePath(kv.first, path));
        }
    }

    std::sort(pathes.begin(), pathes.end());
    pathes.erase(std::unique(pathes.begin(), pathes.end()), pathes.end());
    std::reverse(pathes.begin(), pathes.end());
    if (pathes.size() > _K)
    {
        pathes.resize(_K);
    }

    return pathes;
}

bool KShortestPath::_fast_pruning(unsigned s, unsigned v, unsigned distance) const
{
    (void)distance;
    const std::map<unsigned, kLabel >& map = _dist_label.at(s);

    auto ite = map.find(v);
    if (ite == map.end())
    {
        return false;
    }
    const kLabel& label = ite->second;
    size_t size = 0;
    for (const auto& kv : label)
    {
        size += kv.second.size();
    }
    return size >= _K;
}

bool KShortestPath::_pruning(unsigned s, unsigned v, unsigned distance) const
{
    DistancePathes pathes = _query(s, v);
    if (pathes.size() < _K)
    {
        return false;
    }
    return distance >= pathes.back().distance;
}

size_t KShortestPath::_get_path_memory(const kPath& path) const
{
    return path.size() * sizeof(unsigned);
}

size_t KShortestPath::_get_pathes_memory(const kPathes& pathes) const
{
    size_t size = 0;
    for (const kPath& path : pathes)
    {
        size += _get_path_memory(path);
    }
    return size;
}

size_t KShortestPath::_get_label_memory(const kLabel& label) const
{
    size_t size = 0;
    for (const auto& kv : label)
    {
        size += sizeof(unsigned);              // key
        size += _get_pathes_memory(kv.second); // value
    }
    return size;
}

size_t KShortestPath::get_memory_usage() const
{
    size_t size = 0;
    for (const kLabel& label : _loop_label)
    {
        size += _get_label_memory(label);
    }
    for (const std::map<unsigned, kLabel >& map : _dist_label)
    {
        for (const auto& kv : map)
        {
            size += sizeof(unsigned);
            size += _get_label_memory(kv.second);
        }
    }
    return size;
}

rapidjson::Value _to_json(const Point& p, rapidjson::Document& document)
{
    rapidjson::Value v(rapidjson::kArrayType);
    v.Reserve(2, document.GetAllocator());
    v.PushBack(p.x, document.GetAllocator());
    v.PushBack(p.y, document.GetAllocator());

    return v;
}

rapidjson::Value _to_json(const std::vector<Point>& points, rapidjson::Document& document)
{
    rapidjson::Value v(rapidjson::kArrayType);
    v.Reserve(points.size(), document.GetAllocator());

    for (const Point& p : points)
    {
        v.PushBack(_to_json(p, document), document.GetAllocator());
    }

    return v;
}

rapidjson::Value _to_json(const std::map<Point, unsigned>& map, rapidjson::Document& document)
{
    rapidjson::Value v(rapidjson::kArrayType);
    v.Reserve(map.size() * 2, document.GetAllocator());

    for (const auto& kv : map)
    {
        v.PushBack(_to_json(kv.first, document), document.GetAllocator());
        v.PushBack(kv.second, document.GetAllocator());
    }

    return v;
}

rapidjson::Value _to_json(const kPath& path, rapidjson::Document& document)
{
    rapidjson::Value v(rapidjson::kArrayType);
    v.Reserve(path.size(), document.GetAllocator());

    for (unsigned p : path)
    {
        v.PushBack(p, document.GetAllocator());
    }

    return v;
}

rapidjson::Value _to_json(const kPathes& pathes, rapidjson::Document& document)
{
    rapidjson::Value v(rapidjson::kArrayType);
    v.Reserve(pathes.size(), document.GetAllocator());

    for (const kPath& p : pathes)
    {
        v.PushBack(_to_json(p, document), document.GetAllocator());
    }

    return v;
}

rapidjson::Value _to_json(const kLabel& label, rapidjson::Document& document)
{
    rapidjson::Value v(rapidjson::kArrayType);
    v.Reserve(label.size() * 2, document.GetAllocator());

    for (const auto& kv : label)
    {
        v.PushBack(kv.first, document.GetAllocator());
        v.PushBack(_to_json(kv.second, document), document.GetAllocator());
    }

    return v;
}

rapidjson::Value _to_json(const std::vector<kLabel>& labels, rapidjson::Document& document)
{
    rapidjson::Value v(rapidjson::kArrayType);
    v.Reserve(labels.size(), document.GetAllocator());

    for (const kLabel& label : labels)
    {
        v.PushBack(_to_json(label, document), document.GetAllocator());
    }

    return v;
}

rapidjson::Value _to_json(const std::map<unsigned, kLabel>& map, rapidjson::Document& document)
{
    rapidjson::Value v(rapidjson::kArrayType);
    v.Reserve(map.size() * 2, document.GetAllocator());

    for (const auto& kv : map)
    {
        v.PushBack(kv.first, document.GetAllocator());
        v.PushBack(_to_json(kv.second, document), document.GetAllocator());
    }

    return v;
}

rapidjson::Value _to_json(const std::vector<std::map<unsigned, kLabel> >& maps, rapidjson::Document& document)
{
    rapidjson::Value v(rapidjson::kArrayType);
    v.Reserve(maps.size(), document.GetAllocator());

    for (const auto& map : maps)
    {
        v.PushBack(_to_json(map, document), document.GetAllocator());
    }

    return v;
}

Point _json_to_point(const rapidjson::Value& v)
{
    assert(v.IsArray());

    Point p(v[0].GetInt(), v[1].GetInt());

    return p;
}

std::vector<Point> _json_to_points(const rapidjson::Value& v)
{
    std::vector<Point> points;
    points.reserve(v.GetArray().Size());

    for (const auto& p : v.GetArray())
    {
        points.push_back(_json_to_point(p));
    }

    return points;
}

std::map<Point, unsigned> _json_to_point_map(const rapidjson::Value& v)
{
    std::map<Point, unsigned> map;

    for (rapidjson::Document::ConstValueIterator ite = v.Begin(); ite != v.End(); ++ite)
    {
        const rapidjson::Value& key = *ite;
        const rapidjson::Value& value = *(++ite);
        map[_json_to_point(key)] = value.GetUint();
    }

    return map;
}

kPath _json_to_path(const rapidjson::Value& v)
{
    kPath path;
    path.reserve(v.GetArray().Size());

    for (const auto& p : v.GetArray())
    {
        path.push_back(p.GetUint());
    }

    return path;
}

kPathes _json_to_pathes(const rapidjson::Value& v)
{
    kPathes pathes;
    pathes.reserve(v.GetArray().Size());

    for (const auto& p : v.GetArray())
    {
        pathes.push_back(_json_to_path(p));
    }

    return pathes;
}

kLabel _json_to_label(const rapidjson::Value& v)
{
    kLabel label;

    for (rapidjson::Document::ConstValueIterator ite = v.Begin(); ite != v.End(); ++ite)
    {
        const rapidjson::Value& key = *ite;
        const rapidjson::Value& value = *(++ite);
        label[key.GetUint()] = _json_to_pathes(value);
    }

    return label;
}

std::vector<kLabel> _json_to_loop_label(const rapidjson::Value& v)
{
    std::vector<kLabel> loop_label;
    loop_label.reserve(v.GetArray().Size());

    for (const auto& p : v.GetArray())
    {
        loop_label.push_back(_json_to_label(p));
    }

    return loop_label;
}

std::map<unsigned, kLabel> _json_to_label_map(const rapidjson::Value& v)
{
    std::map<unsigned, kLabel> label_map;

    for (rapidjson::Document::ConstValueIterator ite = v.Begin(); ite != v.End(); ++ite)
    {
        const rapidjson::Value& key = *ite;
        const rapidjson::Value& value = *(++ite);
        label_map[key.GetUint()] = _json_to_label(value);
    }

    return label_map;
}

std::vector<std::map<unsigned, kLabel> > _json_to_dist_label(const rapidjson::Value& v)
{
    std::vector<std::map<unsigned, kLabel> > dist_label;
    dist_label.reserve(v.GetArray().Size());

    for (const auto& p : v.GetArray())
    {
        dist_label.push_back(_json_to_label_map(p));
    }

    return dist_label;
}

void KShortestPath::read_index(const char* file)
{
    FILE* fp = fopen(file, "r");

    char readBuffer[65536];
    rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
    rapidjson::Document d;
    if (d.ParseStream(is).HasParseError())
    {
        printf("%s\n", rapidjson::GetParseError_En(d.GetParseError()));
        return;
    }

    assert(d.HasMember("K"));
    assert(d.HasMember("V"));
    assert(d.HasMember("points"));
    assert(d.HasMember("point_map"));
    assert(d.HasMember("loop_label"));
    assert(d.HasMember("dist_label"));
    _K = d.FindMember("K")->value.GetInt();
    _V = d.FindMember("V")->value.GetInt();
    _points = _json_to_points(d.FindMember("points")->value);
    _point_map = _json_to_point_map(d.FindMember("point_map")->value);
    _loop_label = _json_to_loop_label(d.FindMember("loop_label")->value);
    _dist_label = _json_to_dist_label(d.FindMember("dist_label")->value);

    fclose(fp);
}

void KShortestPath::write_index(const char* file) const
{
    rapidjson::Document d(rapidjson::kObjectType);

    d.AddMember("K", _K, d.GetAllocator());
    d.AddMember("V", _V, d.GetAllocator());
    d.AddMember("points", _to_json(_points, d), d.GetAllocator());
    d.AddMember("point_map", _to_json(_point_map, d), d.GetAllocator());
    d.AddMember("loop_label", _to_json(_loop_label, d), d.GetAllocator());
    d.AddMember("dist_label", _to_json(_dist_label, d), d.GetAllocator());

    FILE* fp = fopen(file, "w");
    char writeBuffer[65536];
    rapidjson::FileWriteStream os(fp, writeBuffer, sizeof(writeBuffer));
    rapidjson::Writer<rapidjson::FileWriteStream> writer(os);

    d.Accept(writer);

    fclose(fp);
}

void KShortestPath::add_indexing_time(double t) {
    _indexing_time += t;
}

void KShortestPath::add_query_time(double t) {
    _query_time += t;
}


void KShortestPath::read_type(const char* map_file)
{
    std::string s(map_file);
    if (s.empty())
    {
        return;
    }
    _read_type(map_file);
}

void KShortestPath::_read_type(const char* map_file)
{
    std::ifstream infile(map_file);
    std::string line;
    unsigned type = 4;
    int length;
    int width;
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

        _insert_map(Point(x, y), type);
    }
}

void KShortestPath::_insert_map(const Point& p, unsigned type)
{
    _type_map.insert(std::make_pair(p, type));
}
