#ifndef SearchState_H
#define SearchState_H
#include <vector>
#include "LocoState.h"
#include "Tools.h"

class SearchState;

struct Neighbor
{
  SearchState* st;
  float length;
  float cost;
};

class SearchState : public LocoState
{

public:
    SearchState();
    SearchState(LocoState);
	~SearchState();

	std::vector<Neighbor> get_virtual_neighbors(bool);
	std::vector<Neighbor> get_complete_neighbors();

	float getCost();
	void setCost(float);
	float getLength();
	void setLength(float);
	SearchState* getPredecessor();
    void setPredecessor(SearchState*);
    //void check_tag(int);

    void insert_virtual_neighbor(SearchState*, float, std::pair<float, float>);
    void insert_complete_neighbor(SearchState*, float, float);

	float getLambda(float, bool);

    bool goalCheck(const SearchState*, float);

    void prune();
    void unprune();
    bool if_pruned();
    void clock();
    void cunlock();
    void llock();
    void lunlock();
    bool if_locked();

    void put_in_frontier();
    bool if_in_frontier();

    void reset_values();
    void reset_status();

	float getCCost();
	void setCCost(float);
	float getCLength();
	void setCLength(float);
	float getLCost();
	void setLCost(float);
	float getLLength();
	void setLLength(float);
	float getLambdaCost();
	void setLambdaCost(float);
	float getLambdaLength();
	void setLambdaLength(float);

	loco_info getLocoInfo();

	void setToStart();

    //SearchState copySearchState();
	
private:
    float cost;
    float length;
	std::vector<Neighbor> virtual_neighbors_alpha;
	std::vector<Neighbor> virtual_neighbors_beta;
	std::vector<Neighbor> complete_neighbors;
	SearchState* predecessor_ptr;

	float c_cost;
	float c_length;
	float l_cost;
	float l_length;
	float lambda_cost;
	float lambda_length;

	//int rand_tag;
	bool pruned;
	bool clocked;
	bool llocked;
	bool infrontier;


};
#endif // SearchState_H