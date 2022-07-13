#ifndef SearchSpace_H
#define SearchSpace_H
#include <vector>
#include <tuple>
#include <algorithm>
#include "SearchState.h"
#include "VirIndex.h"
#include "PhyIndex.h"

//Shao-Heng's definitions
struct solution
{
  std::vector<SearchState*> path;
  float path_length;
  float path_cost;
};

struct heuristic
{
    float min_length;
    float cost_of_min_length;
    float min_cost;
    float length_of_min_cost;
};

struct larac_result
{
  solution sol;
  float lambda;
  float obj_lower_bound;
};

struct fptas_backtrace
{
    float cost;
    SearchState* predecessor;
    int steplength;
};

class SearchSpace
{

public:
    SearchSpace(int, PhyIndex&);
	~SearchSpace();

	bool build_states(const char*, const char*, std::string&, std::string&, bool = false); //take virtual index, physical index

	//subroutines of build_states
	bool construct_virtual_states(const char*);
	int read_physical_index_and_alpha_beta(std::ifstream&, std::string&);
	bool construct_complete_states(const char*, bool = false);
	bool build_all_neighbors(std::ifstream&, int, bool = false);

	float prune_space(SearchState*, SearchState*, float, bool = false, float = 0);
	std::vector<solution> find_sol(SearchState*, SearchState*, float, float, bool = false, float = 0);

	float decide_interval(float, float, float);

	void prune_locked_and_unvisited();
	void reset_status(bool virtual_flag = false);
	void reset_values(bool virtual_flag = false);
	void clear_map();

	solution heuristic_search(SearchState*, SearchState*, bool = false, bool = false, bool = false, float = 0, float = std::numeric_limits<int>::max(), bool = false);
	larac_result larac_search(SearchState*, SearchState*, float, float, bool = false);

	SearchState* getLeader(SearchState*);
	heuristic getHeuristic(SearchState* state, SearchState* dest, bool alpha_or_beta = true);
	heuristic gen_heuristic(SearchState* state, SearchState* dest, bool alpha_or_beta = true);
	//void gen_all_heuristic();

	solution Larac_Advised_Astar(SearchState* startState, SearchState* goalState, float c, float eps = 0.1, bool filter = false, bool user_study = false);
	solution pruning_and_fptas(SearchState* startState, SearchState* goalState, float c, float eps = 0.1, bool pruning_flag = true, bool lambda_flag = false, bool user_study = false );
	solution cost_first(SearchState* startState, SearchState* goalState, float c, float eps = 0.1, bool pruning_flag = true, bool lambda_flag = false, bool user_study = false );
	PhyPath optimize_solution(solution& sol, SearchState* startState);

	SearchState* getState(loco_info);
	SearchState* getVirState(loco_info);

	std::pair<float, float> get_alpha_beta_DP(float);

private:
	std::map<loco_info, SearchState*> state_map;
	std::map<loco_info, SearchState*> vir_state_map;
	VirIndex v;
	PhyIndex* p_ptr;

	std::map<std::pair<SearchState*, SearchState*>, std::pair<heuristic, heuristic> > heuristic_map;
	std::vector<float> _sh_alpha_labels;
	std::vector<float> _sh_beta_labels;

	bool heuristics;

};
#endif // SearchSpace_H