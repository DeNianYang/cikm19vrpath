#include "SearchState.h"

SearchState::SearchState()
{

}

SearchState::SearchState(LocoState lc):LocoState(lc)
{
    reset_values();
    predecessor_ptr = NULL;
}

/*
SearchState SearchState::copySearchState() {
    // TODO
}
*/

void SearchState::insert_virtual_neighbor(SearchState* st_ptr, float l, std::pair<float,float> c)
{
    Neighbor n_alpha = {st_ptr, l, c.first};
    virtual_neighbors_alpha.push_back(n_alpha);

    Neighbor n_beta = {st_ptr, l, c.second};
    virtual_neighbors_beta.push_back(n_beta);
}

void SearchState::insert_complete_neighbor(SearchState* st_ptr, float l, float c)
{
    Neighbor n = {st_ptr, l, c};
    complete_neighbors.push_back(n);
}


std::vector<Neighbor> SearchState::get_virtual_neighbors(bool alpha_or_beta)
{
    if( !(virtual_neighbors_alpha.empty()) ){
        if ( alpha_or_beta ) { return virtual_neighbors_alpha; }
    	else { return virtual_neighbors_beta; }
    }
    else{
        if ( alpha_or_beta ) { return virtual_neighbors_alpha; }
        else { return virtual_neighbors_beta; }
    	//TODO
    }
}

std::vector<Neighbor> SearchState::get_complete_neighbors()
{
    if( !(complete_neighbors.empty()) ){
    	return complete_neighbors;
    }
    else{
        return complete_neighbors;
    	//TODO
    }
}


bool SearchState::goalCheck(const SearchState* st, float c)
{
    return ((gamma_v_x == st->gamma_v_x) && (gamma_v_y == st->gamma_v_y) && (cost <= c));
}

float SearchState::getCost(){ return cost; }
void SearchState::setCost(float c){ cost = c; }
float SearchState::getLength(){ return length; }
void SearchState::setLength(float l){ length = l; }
SearchState* SearchState::getPredecessor(){ return predecessor_ptr; }
void SearchState::setPredecessor(SearchState* st){ predecessor_ptr = st; }

/*void SearchState::check_tag(int tag){
    if(tag != rand_tag){
        rand_tag = tag;
        cost = std::numeric_limits<float>::max();
        length = std::numeric_limits<float>::max();
        //std::cout << "Reseting values for state " << gamma_v_x << std::endl;
    }
}*/

float SearchState::getLambda(float lambda, bool flag = false){
    //print();
    //std::cout << "getLambda called, cost = " << cost << ", length = " << length << ", lambda = " << lambda << std::endl;
	if (flag) {return cost;}
    else {return length + lambda * cost;}
}

void SearchState::prune(){
    pruned = true;
    #ifdef DEBUG
    std::cout << "Pruning state";
    print();
    #endif //DEBUG
}
void SearchState::unprune(){ pruned = false; }
bool SearchState::if_pruned() { return pruned; }

void SearchState::clock(){
    clocked = true;
    #ifdef DEBUG
    std::cout << "C-Locking state";
    print();
    #endif //DEBUG
}
void SearchState::cunlock(){ clocked = false; }
void SearchState::llock(){
    llocked = true;
    #ifdef DEBUG
    std::cout << "L-Locking state";
    print();
    #endif //DEBUG
}
void SearchState::lunlock(){ llocked = false; }
bool SearchState::if_locked() { return (clocked || llocked); }
bool SearchState::if_in_frontier() { return infrontier; }
void SearchState::put_in_frontier() { infrontier = true; }

void SearchState::reset_values() {
    float infty = std::numeric_limits<float>::max();
    cost = infty;
    length = infty;
    c_cost = infty;
    c_length = infty;
    l_cost = infty;
    l_length = infty;
    lambda_cost = infty;
    lambda_length = infty;
}

void SearchState::reset_status() {
    clocked = false;
    llocked = false;
    infrontier = false;
    pruned = false;
}

float SearchState::getCCost(){ return c_cost; }
void SearchState::setCCost(float c){ c_cost = c; }
float SearchState::getCLength(){ return c_length; }
void SearchState::setCLength(float l){ c_length = l; }
float SearchState::getLCost(){ return l_cost; }
void SearchState::setLCost(float c){ l_cost = c; }
float SearchState::getLLength(){ return l_length; }
void SearchState::setLLength(float l){ l_length = l; }
float SearchState::getLambdaCost(){ return lambda_cost; }
void SearchState::setLambdaCost(float c){ lambda_cost = c; }
float SearchState::getLambdaLength(){ return lambda_length; }
void SearchState::setLambdaLength(float l){ lambda_length = l; }

loco_info SearchState::getLocoInfo()
{
    loco_info lc = {gamma_v_x, gamma_v_y, theta_v, gamma_p_x, gamma_p_y, theta_p};
    return lc;
}

void SearchState::setToStart(){
    setCost(0);
    setLength(0);
    setCCost(0);
    setCLength(0);
    setLCost(0);
    setLLength(0);
    setLambdaCost(0);
    setLambdaLength(0);
}

SearchState::~SearchState()
{
}
