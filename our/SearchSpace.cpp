#include "SearchSpace.h"
#include <queue>
#include <limits>
#include <cmath>
#include <ctime>

//#define PRUNEDEBUG
#define DEBUGG
//#define DEBUG_MCPP

bool operator <(const loco_info& x, const loco_info& y) {
    return std::tie(x.gamma_v_x, x.gamma_v_y, x.theta_v, x.gamma_p_x, x.gamma_p_y, x.theta_p)
    < std::tie(y.gamma_v_x, y.gamma_v_y, y.theta_v, y.gamma_p_x, y.gamma_p_y, y.theta_p);
}

// the top of the priority queue is the greatest element by default, but we want the smallest, so flip the sign
bool operator<(const std::pair<SearchState*, float> st1, const std::pair<SearchState*, float> st2) {
  return st1.second > st2.second;
}

int downscale(float value, float interval){
  return static_cast<int>(std::ceil(value/interval));
}

float upscale(float value, float interval){
  return value * interval;
}

/*void PhyPrint(PhyPath& p, std::vector<Step>& vsteps, Point startpoint){
    for (unsigned i = 0; i <= vsteps.size(); ++i){
      if (i <= vsteps.size() ){
        std::cout << '(' << startpoint.x << ',' << startpoint.y << ")";
      }
      else{
        std::cout << "(  ,  )";
      }
      if (i < vsteps.size() ){
        startpoint.x += vsteps[i].delta_x;
        startpoint.y += vsteps[i].delta_y;
      }
      if (i < p.physical_path.size() ){
        std::cout << " | (" << p.physical_path[i].x << ',' << p.physical_path[i].y << ") | ";
      }
      else{
        std::cout << " | ( , ) | ";
      }
      if (i < p.vr_ops.size() ){
        std::cout << p.vr_ops[i];
      }
      std::cout << std::endl;
    }
    std::cout << "Distance = " << p.distance << ", Cost = " << p.cost << std::endl;
}*/

SearchSpace::SearchSpace(int max_lv, PhyIndex& p):v(max_lv)
{
	heuristics = false;
  p_ptr = &p;
}

SearchSpace::~SearchSpace(){
  clear_map();
}

bool SearchSpace::build_states(const char* vir_graph, const char* phy_graph, std::string& vir_index, std::string& phy_index, bool simplify)
{
  clock_t t = clock();

  std::ifstream pfile(phy_index);
  int table_size = -1;

  bool flag = false;
  flag = construct_virtual_states(vir_graph);
  if (flag) {
    table_size = read_physical_index_and_alpha_beta(pfile, vir_index);
  }

  flag = false;
  if (table_size != -1) { flag = construct_complete_states(phy_graph, simplify); }
  if (flag) { flag = build_all_neighbors(pfile, table_size, simplify); }
  if (flag) {
    std::cout << "Search Space Built in " << float(clock() - t)/CLOCKS_PER_SEC << " seconds.\n";
    t = clock();
    //gen_all_heuristic();
    //std::cout << "Heuristics Built in " << float(clock() - t)/CLOCKS_PER_SEC << " seconds.\n";

    #ifdef DEBUG
    for(unsigned i = 0; i < _sh_alpha_labels.size(); ++i){
      std::cout << "_sh_alpha_labels[" << i << "] = " << _sh_alpha_labels[i] << std::endl;
    }
    for(unsigned i = 0; i < _sh_beta_labels.size(); ++i){
      std::cout << "_sh_beta_labels[" << i << "] = " << _sh_beta_labels[i] << std::endl;
    }
    #endif // DEBUG

    return true;
  }
  return false;

}

bool SearchSpace::construct_virtual_states(const char* vir_graph)
{
  std::cout << "Constructing Virtual States...\n";
  //Read the virtual input graph to get virtual point map
  v.read_graph(vir_graph);
  std::map<Point, int>* v_map = v.get_point_map();

  if ( v_map->size() == 0) { return false; }

  //Put the elements in v_map into vir_state_map
  for ( auto it = v_map->begin(); it != v_map->end(); it++ ){
    if (it->second == 3) { continue; }// is obj
    int x = it->first.x;
    int y = it->first.y;
    //for now, virtual map does not contain angles
    loco_info lc = {x, y, 0, 0, 0, 0};
    SearchState* st = new SearchState(LocoState(x,y,0,0,0,0));
    vir_state_map[lc] = st;
  }

  return true;
}

int SearchSpace::read_physical_index_and_alpha_beta(std::ifstream& pfile, std::string& vir_index) //returns table size = lmax + 1
{
  std::cout << "Reading Physical Index...\n";
  //Read alpha, beta values from physical index

  std::string line;
  std::getline(pfile, line);
  std::istringstream iss(line);
  int lmax = 0;
  if ( !(iss >> lmax) ) {
    std::cout << "Read Physical Index Error: missing l_max" << "\n";
    return -1;
  }
  p_ptr->set_max_lv(lmax);

  std::getline(pfile, line);
  if (line != "ab") {
    std::cout << "Read Physical Index Error: missing a-b values" << "\n";
    return -1;
  }

  int table_size = lmax + 1;

  _sh_alpha_labels.push_back(0);
  _sh_beta_labels.push_back(0); // initialize alpha(0) = 0, beta(0) = 0;

  for (int i = 1; i < table_size; ++i){
    std::getline(pfile, line);
    std::istringstream iss(line);
    float a, b;
    if( !(iss >> a >> b) ) {
      std::cout << "Read Physical Index Error: missing values at l = " << i << "\n";
      return -1;
    }
    _sh_alpha_labels.push_back(a);
    _sh_beta_labels.push_back(b);
  }

  //Read the virtual index, build neighbors
  std::ifstream vfile(vir_index);
  while (std::getline(vfile, line)){
    if (line == "---") {
      std::string xyline;
      int x, y;
      std::getline(vfile, xyline);
      std::istringstream xystream(xyline);
      if (!(xystream >> x >> y)) {
        std::cout << "Read Virtual Index Error: missing point value" << "\n";
        return -1;
      }
      Point start(x,y);
      loco_info lc1 = {x, y, 0, 0, 0, 0};

      std::getline(vfile, xyline);
      std::istringstream xystream2(xyline);
      if (!(xystream2 >> x >> y)) {
        std::cout << "Read Virtual Index Error: missing point value" << "\n";
        return -1;
      }
      Point end(x,y);
      loco_info lc2 = {x, y, 0, 0, 0, 0};

      std::getline(vfile, xyline);
      std::istringstream xystream3(xyline);
      if (!(xystream3 >> x >> y)) {
        std::cout << "Read Virtual Index Error: missing point value" << "\n";
        return -1;
      }

      auto start_it = vir_state_map.find(lc1);
      auto end_it = vir_state_map.find(lc2);
      if ( !(start_it == vir_state_map.end()) && !(end_it == vir_state_map.end()) ){

        float length = std::pow( std::pow(x, 2.0) + std::pow(y, 2.0), 0.5);
        std::pair<float, float> abcost = get_alpha_beta_DP(length);

        if (abcost.first < 1000){
          start_it->second->insert_virtual_neighbor(end_it->second, length, abcost);

          #ifdef DEBUG
          std::cout << "Inserting virtual neighbor ";
          end_it->second->print();
          std::cout << " to ";
          start_it->second->print();
          std::cout << ", length = " << length << ", cost = " << abcost.first << "," << abcost.second << std::endl;
          #endif //DEBUG
        }
      }
      else {
        std::cout << "Read Virtual Index Error: Point(s) not in virtual map" << "\n";
        return -1;
      }
    }
  }

  #ifdef DEBUG
  std::cout << "Virtual map size = " << vir_state_map.size() << std::endl;
  #endif // DEBUG

  return table_size;
}

std::pair<float, float> SearchSpace::get_alpha_beta_DP(float length){
  unsigned ll = std::floor(length);
  unsigned ul = std::ceil(length);

  if (_sh_alpha_labels.size() == 1) {
    std::cout << "Error: get_alpha_beta_DP called with an empty alpha-beta table\n";
    return std::make_pair(0,0);
  }
  else if (_sh_alpha_labels.size() >= ul + 1) { //table is ready for this query
    return std::make_pair(_sh_alpha_labels[ll], _sh_beta_labels[ul]);
  }
  else{ // do not generate edges larger than lmax!
    /*
    int nowmax = _sh_alpha_labels.size() - 1; //get the max column
    for (unsigned i = nowmax + 1; i <= ul; ++i){
      _sh_alpha_labels.push_back(std::numeric_limits<float>::max());
      _sh_beta_labels.push_back(0);

      for (int j = 1; j <= (i * 0.5); ++j){
        if ( (_sh_alpha_labels[j] + _sh_alpha_labels[i-j]) < _sh_alpha_labels[i] ) {
          _sh_alpha_labels[i] = _sh_alpha_labels[j] + _sh_alpha_labels[i-j];
        }
        if ( (_sh_beta_labels[j] + _sh_beta_labels[i-j]) > _sh_beta_labels[i] ) {
          _sh_beta_labels[i] = _sh_beta_labels[j] + _sh_beta_labels[i-j];
        }
      }
    }
    return std::make_pair(_sh_alpha_labels[ll], _sh_beta_labels[ul]);
    */
    return std::make_pair(9999,9999);
  }
}

bool SearchSpace::construct_complete_states(const char* phy_graph, bool simplify)
{
  std::cout << "Constructing Complete States...\n";
  p_ptr->read_graph(phy_graph); //Read the physical input graph to get physical point map
  std::map<Point, int>* p_map = p_ptr->get_point_map();

  if ( p_map->size() == 0 ) {
    std::cout << "Read Physical Input Error\n";
    return false;
  }

  if ( vir_state_map.size() == 0 ) {
    std::cout << "Load Virtual State Map Error\n";
    return false;
  }

  //for states in vir_state_map, pair it with all physical positions and throw into state_map
  for (auto it1 = vir_state_map.begin(); it1 != vir_state_map.end(); ++it1){
    int vx = it1->first.gamma_v_x;
    int vy = it1->first.gamma_v_y;

    for (auto it2 = p_map->begin(); it2 != p_map->end(); ++it2){
      int px = it2->first.x;
      int py = it2->first.y;

      //insert one physical state
      loco_info lc = {vx, vy, 0, px, py, 0};
      SearchState* st = new SearchState(LocoState(vx, vy, 0, px, py, 0));
      state_map[lc] = st;

      if (!simplify){ //insert other seven physical states.
          for(int i = 45; i < 360; i += 45){
          loco_info lc = {vx, vy, 0, px, py, i};
          SearchState* st = new SearchState(LocoState(vx, vy, 0, px, py, i));
          state_map[lc] = st;
        }
      }
    }
  }

  #ifdef DEBUG
  std::cout << "State map size = " << state_map.size() << std::endl;
  #endif // DEBUG

  return true;
}

bool SearchSpace::build_all_neighbors(std::ifstream& pfile, int table_size, bool simplify)
{
  std::cout << "Building All Neighbors...\n";
  //Build ShaoHeng's physical index from file
  std::string line;
  std::getline(pfile, line);
  if (line != "Index") {
    std::cout << "Build Neighbors Error: no SH's index detected\n";
    return false;
    }
  int p1x, p1y, p1theta, p2x, p2y;

  while(std::getline(pfile, line)) {
    std::istringstream iss(line);
    if ( !(iss >> p1x >> p1y >> p1theta >> p2x >> p2y ) ) {
      std::cout << "Build Neighbors Error: Incomplete State\n";
      return false;
    }

    std::getline(pfile, line);
    std::string line2 = line;
    std::istringstream indexss(line);
    float indexed_cost;

    std::vector<float> cost_table;

    for(int i = 0; i < table_size; ++i){
      if ( !( indexss >> indexed_cost) ) {
        std::cout << "Table size: " << table_size << std::endl;
        std::cout << line2;
        std::cout << "Build Neighbors Error: Incomplete Index\n";
        return false;
      }
      cost_table.push_back(indexed_cost);
    }

    for (auto it1 = vir_state_map.begin(); it1 != vir_state_map.end(); ++it1){ //start virtual state
      loco_info lcv1 = it1->first;
      std::vector<Neighbor> neighbors = it1->second->get_virtual_neighbors(true); //does not matter alpha or beta
      for (auto it2 = neighbors.begin(); it2 != neighbors.end(); ++it2){
        unsigned l = std::ceil(it2->length);

        if ( l >= cost_table.size() ) { continue; } //path is too long
        if ( cost_table[l] == -1 ) { continue; } // no path
        SearchState* st2 = it2->st;
        loco_info lcp1 = {lcv1.gamma_v_x, lcv1.gamma_v_y, lcv1.theta_v, p1x, p1y, p1theta};

        // try insert neighbor with angle = 0
        loco_info lcp2 = {st2->getGammaVX(), st2->getGammaVY(), st2->getThetaV(), p2x, p2y, 0};
        auto p1 = state_map.find(lcp1);
        auto p2 = state_map.find(lcp2);
        if ( !(p1 == state_map.end() ) && !(p2 == state_map.end() ) ){
          p1->second->insert_complete_neighbor(p2->second, it2->length, cost_table[l]);
        }

        if (!simplify) { // try insert the other seven neighbors
          for (int angle = 45; angle < 360; angle += 45){
            loco_info lcp2 = {st2->getGammaVX(), st2->getGammaVY(), st2->getThetaV(), p2x, p2y, angle};
            auto p1 = state_map.find(lcp1);
            auto p2 = state_map.find(lcp2);
            if ( !(p1 == state_map.end() ) && !(p2 == state_map.end() ) ){ p1->second->insert_complete_neighbor(p2->second, it2->length, cost_table[l]); }
          }
        }
      }
    }
  }
  return true;
}

SearchState* SearchSpace::getLeader(SearchState* st)
{
  loco_info info = {st->getGammaVX(), st->getGammaVY(), st->getThetaV(), 0, 0, 0};
  auto iter = vir_state_map.find(info);
  if (iter!= vir_state_map.end()){
    return iter->second;
  }
  else {

    #ifdef DEBUG
      std::cout << "Error: Cannot find the virtual leader for search state ";
      st->print();
    #endif //DEBUG

    return NULL;
  }
}

heuristic SearchSpace::getHeuristic(SearchState* currentState, SearchState* goalState, bool alpha_or_beta)
{
  //need to use their virtual leaders to get values
  auto leader_currentState = getLeader(currentState);
  auto leader_goalState = getLeader(goalState);
  if ((leader_currentState != NULL) && (leader_goalState != NULL)){
    auto iter = heuristic_map.find(std::make_pair(leader_currentState, leader_goalState));
    if (iter != heuristic_map.end()){
      return (alpha_or_beta? iter->second.first : iter->second.second );
    }
    else {
      return gen_heuristic(leader_currentState, leader_goalState, alpha_or_beta);
    }
  }
}

solution SearchSpace::heuristic_search(SearchState* startState, SearchState* goalState,
  bool heuristic_flag, bool virtual_flag, bool cost_flag, float lambda, float c, bool alpha_or_beta)
{

  if (c != 0) { lambda += 0.01; } // small fix: so that the ub path and the lb path will not be equal
  else { c = std::numeric_limits<float>::max(); } // so no constraint
  reset_values(virtual_flag); // reset the cost / length values

  std::priority_queue<std::pair<SearchState*, float> > frontier;
  std::vector<SearchState*> opt_path;
  startState->setToStart();
  frontier.push(std::make_pair(startState,0));

  while (!frontier.empty()) {
    SearchState* currentState = frontier.top().first; // .top() doesn't actually remove the node
    frontier.pop();
    if ( currentState->if_pruned() ) {continue;} // if this state is pruned, no need to check it again
    currentState->prune(); // use the prune tag to mark as visited

    #ifdef DEBUG
      std::cout << "currentState is ";
      currentState->print();
      std::cout << "Cost = " << currentState->getCost() << std::endl;
    #endif //DEBUG

    if (currentState->goalCheck(goalState, c)) {

      #ifdef DEBUG
        std::cout << "Goal Checked! Cost = " << currentState->getCost() << " , Length = " << currentState->getLength() <<std::endl;
      #endif //DEBUG

      while(currentState != startState){
        opt_path.push_back(currentState);
        currentState = currentState->getPredecessor();
      }

      opt_path.push_back(startState);
      break;
    }

    std::vector<Neighbor> neighbors;

    if(virtual_flag){
      neighbors = currentState->get_virtual_neighbors(alpha_or_beta);
    }
    else{
      neighbors = currentState->get_complete_neighbors();
    }

    while (!neighbors.empty()) {
    	Neighbor current_neighbor = neighbors.back();
    	neighbors.pop_back();
      if ( current_neighbor.st->if_pruned() ) { continue; }

    	float tempcost = currentState->getCost() + current_neighbor.cost;
    	float templength = currentState->getLength() + current_neighbor.length;
    	float templambda;
    	if(cost_flag) {templambda = tempcost;}
    	else {templambda = templength + lambda * tempcost;}

      #ifdef DEBUG
        std::cout << "tempcost = " << tempcost << std::endl;
        std::cout << "templength = " << templength << std::endl;
        std::cout << "templambda = " << templambda << std::endl;
      #endif //DEBUG
      
    	if (templambda < current_neighbor.st->getLambda(lambda, cost_flag)){ //relax
    		current_neighbor.st->setCost(tempcost);
    		current_neighbor.st->setLength(templength);
        current_neighbor.st->setPredecessor(currentState);

        #ifdef DEBUG
        std::cout << "Relaxing state " << current_neighbor.st->getGammaVX() << " : cost = " << tempcost << ", length = " << templength << std::endl;
        #endif //DEBUG
    	}

    	float value = current_neighbor.st->getLambda(lambda, cost_flag);
    	if (heuristic_flag){
        heuristic heu = getHeuristic(current_neighbor.st, goalState, true); //use alpha heuristic
    		value += (heu.min_length + lambda * heu.min_cost);
    	}

        frontier.push(std::make_pair(current_neighbor.st,value));

        #ifdef DEBUG
        std::cout << "Pushing state " << current_neighbor.st->getGammaVX() << " , value = " << value << std::endl;
        #endif //DEBUG
    }
  }

  //find a solution
  if(opt_path.empty()){
    float infty = std::numeric_limits<float>::max();
    solution nosol = {opt_path, infty, infty};
    return nosol;
  }
  else{
    solution sol = {opt_path, opt_path.front()->getLength(), opt_path.front()->getCost()};
    //std::cout << "returning solution with cost = " << opt_path.front()->getCost() << std::endl;
    //std::cout << "This solution path ends at";
    //opt_path.front()->print();
    return sol;
  }
}

larac_result SearchSpace::larac_search(SearchState* startState, SearchState* goalState, float c, float eps, bool alpha_or_beta)
{

  float ratio;
  solution sol_c, sol_l, sol_temp, sol;

  reset_status(true); //is used in virtual map, so virtual_flag = true
  sol_l = heuristic_search(
    startState,
    goalState,
    false, //heuristic_flag
    true, //virtual_flag
    false, //cost_flag
    0, // lambda = 0
    0, // no cost constraint
    alpha_or_beta // alpha_or_beta
    );
  if (sol_l.path_cost <= c){

    #ifdef DEBUG
    std::cout << "Greedy Ready! Lambda = -1" << std::endl;
    #endif //DEBUG

    larac_result result = {sol_l, -1, 0};
    return result;
  }

  reset_status(true); //is used in virtual map, so virtual_flag = true
  sol_c = heuristic_search(
    startState,
    goalState,
    false, //heuristic_flag
    true, //virtual_flag
    true, //cost_flag
    0, // lambda = 0
    0, // no cost constraint
    alpha_or_beta // alpha_or_beta
    );
  if (sol_c.path_cost > c){

    #ifdef DEBUG
    std::cout << "Infeasible!" << std::endl;
    #endif //DEBUG

    larac_result result = {sol_c, -2, std::numeric_limits<float>::max() };
    return result;
  }

  while(true){
    ratio = (sol_l.path_length - sol_c.path_length)/(sol_c.path_cost - sol_l.path_cost);
    reset_status(true); //is used in virtual map, so virtual_flag = true
    sol_temp = heuristic_search(
      startState,
      goalState,
      false, //heuristic_flag
      true, //virtual_flag
      false, //cost_flag
      ratio, //lambda
      0, // no cost constraint
      alpha_or_beta // alpha_or_beta
      );

    if(abs(sol_l.path_cost - sol_c.path_cost)/sol_l.path_cost <= eps){
      larac_result result = {sol_c, ratio, sol_l.path_length };
      return result;
    }
    else if(abs((sol_temp.path_length + ratio * sol_temp.path_cost) == (sol_l.path_length + ratio * sol_l.path_cost))){
      larac_result result = {sol_c, ratio, sol_l.path_length };
      return result;
    }
    else if(sol_temp.path_cost <= c){
      sol_c = sol_temp;
    }
    else{
      sol_l = sol_temp;
    }
  }
}

heuristic SearchSpace::gen_heuristic(SearchState* s, SearchState* t, bool alpha_or_beta)
{
  heuristic ha = {0,0,0,0};
  heuristic hb = ha;
  if (s == t){ //self-to-self
    heuristic_map[std::make_pair(s,s)] = std::make_pair(ha,ha);
    return ha;
  }

  else{
    #ifdef DEBUG
      std::cout << "Preparing to generate heuristics between";
      s->print();
      std::cout << " and ";
      t->print();
      std::cout << std::endl;
    #endif //DEBUG

    //generate alpha heuristics
    reset_status(true); //is used in virtual map, so virtual_flag = true
    solution sol_a_l = heuristic_search(
      s,
      t,
      false, //heuristic_flag
      true, //virtual_flag
      false, //cost_flag
      0, // lambda = 0
      0, // no cost constraint
      true // heuristics take alpha
      );
    reset_status(true); //is used in virtual map, so virtual_flag = true
    solution sol_a_c = heuristic_search(
      s,
      t,
      false, //heuristic_flag
      true, //virtual_flag
      true, //cost_flag
      0, // lambda = 0
      0, // no cost constraint
      true // heuristics take alpha
      );

    //generate beta heuristics
    reset_status(true); //is used in virtual map, so virtual_flag = true
    solution sol_b_l = heuristic_search(
      s,
      t,
      false, //heuristic_flag
      true, //virtual_flag
      false, //cost_flag
      0, // lambda = 0
      0, // no cost constraint
      false // heuristics take beta
      );
    reset_status(true); //is used in virtual map, so virtual_flag = true
    solution sol_b_c = heuristic_search(
      s,
      t,
      false, //heuristic_flag
      true, //virtual_flag
      true, //cost_flag
      0, // lambda = 0
      0, // no cost constraint
      false // heuristics take beta
      );

    ha = {sol_a_l.path_length, sol_a_l.path_cost, sol_a_c.path_cost, sol_a_c.path_length};
    hb = {sol_b_l.path_length, sol_b_l.path_cost, sol_b_c.path_cost, sol_b_c.path_length};
    heuristic_map[std::make_pair(s,t)] = std::make_pair(ha,hb);

    #ifdef DEBUG
      std::cout << "Alpha: Length = " << sol_a_l.path_length << ", Cost = " << sol_a_c.path_cost << std::endl;
      std::cout << "Beta: Length = " << sol_b_l.path_length << ", Cost = " << sol_b_c.path_cost << std::endl;
    #endif //DEBUG
  }
  reset_status(true); //is used in virtual map, so virtual_flag = true

  //return the required version
  return ( alpha_or_beta ? ha : hb );
}

float SearchSpace::prune_space(SearchState* startState, SearchState* goalState, float c, bool lambda_flag, float lambda)
{
	#ifdef DEBUG
  	if(!heuristics){
  		std::cout << "WARNING: heuristics not generated before pruning" << std::endl;
  	}
	#endif //DEBUG

  float eps = 0.1;

  startState->setToStart();

  std::priority_queue<SearchState*> frontier;

  startState->setToStart();
  frontier.push(startState);
  startState->put_in_frontier();

  float best = std::numeric_limits<float>::max();

  while (!frontier.empty()) {
    SearchState* currentState = frontier.top(); // .top() doesn't actually remove the node

    #ifdef PRUNEDEBUG
    std::cout << "Checking state";
    currentState->print();
    #endif //DEBUG

    frontier.pop();

    // if the top is a locked state, then the pruning ends
    if (currentState->if_locked()){
      prune_locked_and_unvisited();
      return best;
    }

    // if this is a goal state, update the best solution
    if (currentState->goalCheck(goalState, std::numeric_limits<float>::max())) { //check Setting C,L,Lambda
      #ifdef PRUNEDEBUG
      std::cout << "Goal State! Best is updated from " << best;
      #endif
      if ( (currentState->getCCost() < c ) && (currentState->getCLength() < best) ) { best = currentState->getCLength(); }
      if ( (currentState->getLCost() < c ) && (currentState->getLLength() < best) ) { best = currentState->getLLength(); }
      if ( (lambda_flag) && (currentState->getLambdaCost() < c) && (currentState->getLambdaLength() < best) ) { best = currentState->getLambdaLength(); }
      #ifdef PRUNEDEBUG
      std::cout << " to" << best << std::endl;
      #endif
      continue; //no need to see the neighbors
    }

    // check if this state can be pruned. Preparing values
    heuristic heu_start_alpha = getHeuristic(startState, currentState, true); //use alpha heuristic
    heuristic heu_goal_alpha = getHeuristic(currentState, goalState, true); //use alpha heuristic
    heuristic heu_start_beta = getHeuristic(startState, currentState, false); //use beta bound
    heuristic heu_goal_beta = getHeuristic(currentState, goalState, false); //use beta bound

    float length_from_start_alpha = heu_start_alpha.min_length;
    float cost_from_start_alpha = heu_start_alpha.min_cost;
    float cost_to_goal_alpha = heu_goal_alpha.min_cost;
    float cost_of_length_to_goal_alpha = heu_goal_alpha.cost_of_min_length;
    float length_to_goal_alpha = heu_goal_alpha.min_length;

    //float length_from_start_beta = heu_start_beta.min_length;
    //float cost_from_start_beta = heu_start_beta.min_cost;
    //float cost_to_goal_beta = heu_goal_beta.min_cost;
    float cost_of_length_to_goal_beta = heu_goal_beta.cost_of_min_length;
    float length_to_goal_beta = heu_goal_beta.min_length;

    #ifdef PRUNEDEBUG
    std::cout << "Length from start = " << length_from_start_alpha << ", Cost from start = " << cost_from_start_alpha << std::endl;
    std::cout << "Length to goal = " << length_to_goal_alpha << ", Cost to goal = " << cost_to_goal_alpha << std::endl;
    #endif //DEBUG

    // Infeasible Pruning
    if ( (cost_from_start_alpha + cost_to_goal_alpha) > c ){
      currentState->prune();
      continue; //no need to see the neighbors
    }

    // Suboptimal Pruning; Relaxed a bit to avoid computing errors
    if ( (length_from_start_alpha + length_to_goal_alpha) > best * (1 + eps) ){
      currentState->prune();
      continue; //no need to see the neighbors
    }

    // State Locking - Cost
    if ( (currentState->getCCost() + cost_to_goal_alpha) > c ){
      currentState->clock(); // no need to push it to the frontier; it will be pushed if unlocked
      continue; //no need to see the neighbors, for now
    }

    // State Locking - Length; Relaxed a bit to avoid computing errors
    if ( (currentState->getLLength() + length_to_goal_alpha) > best * (1 + eps) ){
      currentState->llock(); // no need to push it to the frontier; it will be pushed if unlocked
      continue; //no need to see the neighbors, for now
    }

    // Best Value Updating: Setting C,L,Lambda
    #ifdef PRUNEDEBUG
    std::cout << "Now CLength is " << currentState->getCLength() << ", length_to_goal_beta is " << length_to_goal_beta << ", best is " << best << std::endl;
    std::cout << "Now LLength is " << currentState->getLLength() << ", length_to_goal_beta is " << length_to_goal_beta << ", best is " << best << std::endl;
    #endif //DEBUG

    float temp = currentState->getCLength() + length_to_goal_beta;
    if (( (currentState->getCCost() + cost_of_length_to_goal_beta) < c) && (temp < best)) {
      best = temp;
      #ifdef PRUNEDEBUG
      std::cout << "Updating the objective to be " << temp << std::endl;
      #endif //DEBUG
    }

    temp = currentState->getLLength() + length_to_goal_beta;
    if (( (currentState->getLCost() + cost_of_length_to_goal_beta) < c) && (temp < best)) {
      best = temp;

      #ifdef PRUNEDEBUG
      std::cout << "Updating the objective to be " << temp << std::endl;
      #endif //DEBUG
    }

    if (lambda_flag) {
      temp = currentState->getLambdaLength() + length_to_goal_beta;
      if (( (currentState->getLambdaCost() + cost_of_length_to_goal_beta) < c) && (temp < best)) {
        best = temp;

        #ifdef PRUNEDEBUG
        std::cout << "Updating the objective to be " << temp << std::endl;
        #endif //DEBUG
      }
    }

    //the neighbor states are checked; needs to check complete neighbors
    std::vector<Neighbor> neighbors;
    neighbors = currentState->get_complete_neighbors();

    while (!neighbors.empty()) {
      Neighbor current_neighbor = neighbors.back();
      neighbors.pop_back();

      float edgecost = current_neighbor.cost;
      float edgelength = current_neighbor.length;

      // Setting C
      if ( (currentState->getCCost() + edgecost) < current_neighbor.st->getCCost() ){
        current_neighbor.st->setCCost(currentState->getCCost() + edgecost);
        current_neighbor.st->setCLength(currentState->getCLength() + edgelength);

        #ifdef PRUNEDEBUG
        std::cout << "Setting state ";
        current_neighbor.st->print();
        std::cout << " CCost = " << currentState->getCCost() + edgecost;
        std::cout << " , CLength = " << currentState->getCLength() + edgelength << std::endl;
        #endif //DEBUG

        current_neighbor.st->cunlock(); //unlock
      }

      // Setting L
      if ( (currentState->getLLength() + edgelength) < current_neighbor.st->getLLength() ){
        current_neighbor.st->setLCost(currentState->getLCost() + edgecost);
        current_neighbor.st->setLLength(currentState->getLLength() + edgelength);

        #ifdef PRUNEDEBUG
        std::cout << "Setting state ";
        current_neighbor.st->print();
        std::cout << " LCost = " << currentState->getLCost() + edgecost;
        std::cout << " , LLength = " << currentState->getLLength() + edgelength << std::endl;
        #endif //DEBUG

        current_neighbor.st->lunlock(); //unlock
      }

      // Setting Lambda
      if ( lambda_flag ) {
        if ( (currentState->getLambdaCost() + edgecost + lambda * (currentState->getLLength() + edgelength))
          < (current_neighbor.st->getLambdaCost() + lambda * current_neighbor.st->getLLength()) ){
          current_neighbor.st->setLambdaCost(currentState->getLambdaCost() + edgecost);
          current_neighbor.st->setLambdaLength(currentState->getLambdaLength() + edgelength);
        }
      }

      //push to frontier if the neighbor is not locked, but is not there
      if ( ( !current_neighbor.st->if_locked() ) && ( !current_neighbor.st->if_in_frontier() ) ){
        frontier.push(current_neighbor.st); //all states are created equal
        current_neighbor.st->put_in_frontier();
      }
    }
  }
  prune_locked_and_unvisited();
  return best;
}

void SearchSpace::prune_locked_and_unvisited()
{
  //prune all the locked states
  for( auto it = state_map.begin(); it != state_map.end(); it++){
    if (it->second->if_locked()) { it->second->prune(); }
    if ( !(it->second->if_in_frontier() ) ) { it->second->prune(); }
  }
}

std::vector<solution> SearchSpace::find_sol(SearchState* startState, SearchState* goalState, float c, float interval, bool do_all, float ub)
{
  clock_t t = clock();

  reset_values(); // clear the values
  startState->setToStart();

  std::map< SearchState*, std::vector< fptas_backtrace > > matrix;
  float bestcost = 0;

  for( auto it = state_map.begin(); it != state_map.end(); it++){
    if (!it->second->if_pruned()) { //only those unpruned states go into the matrix
      std::vector< fptas_backtrace > temp;

      bestcost = std::numeric_limits<float>::max();
      if (it->second == startState){ bestcost = 0; } //initialize the base case c=0

      SearchState* nptr = NULL;
      fptas_backtrace fb = {bestcost, nptr, 0};
      temp.push_back(fb);
      matrix[it->second] = temp; // create a column
    }
  }

  //heuristic hbound = getHeuristic(startState, goalState);
  //float cap = hbound.length_of_min_cost;
  int trivial_bound = vir_state_map.size() * (std::max(v.width, v.length));
  if(do_all || (ub==0) || (ub>trivial_bound)){
    ub = trivial_bound;
  }

  std::cout << "find_sol: upper bound set at " << ub << std::endl;

  for(unsigned i = 1; i < ub; i++){ // no end loop; needs to break if find solution

    /*
    if ((clock() - t)/CLOCKS_PER_SEC > 1000){
      throw 99;
      break;
    }
    */

    for( auto it = matrix.begin(); it != matrix.end(); it++){

      bestcost = std::numeric_limits<float>::max();
      SearchState* best_predecessor = NULL;
      int best_steplength = 0;
      SearchState* currentState = it->first;
      std::vector<Neighbor> neighbors = currentState->get_complete_neighbors();

      while (!neighbors.empty()) {
        Neighbor current_neighbor = neighbors.back();
        neighbors.pop_back();

        if (do_all){
          //Point s(currentState->getGammaPX(), currentState->getGammaPY());
          //Point t(current_neighbor.st->getGammaPX(), current_neighbor.st->getGammaPY());
          float vdx = current_neighbor.st->getGammaVX() - currentState->getGammaVX();
          float vdy = current_neighbor.st->getGammaVY() - currentState->getGammaVY();
          float pdx = current_neighbor.st->getGammaPX() - currentState->getGammaPX();
          float pdy = current_neighbor.st->getGammaPY() - currentState->getGammaPY();
          float vl = std::pow((std::pow(vdx, 2) + std::pow(vdy,2)), 0.5);
          float va = std::atan2(vdx, vdy);
          float pl = std::pow((std::pow(pdx, 2) + std::pow(pdy,2)), 0.5);
          float pa = std::atan2(pdx, pdy);
          p_ptr->RW_cost(currentState->getThetaV(), currentState->getThetaP(), vl, va, pl, pa, true);
        }

        if ( !current_neighbor.st->if_pruned() ){
          auto neighbor_it = matrix.find(current_neighbor.st);
          int length = downscale(current_neighbor.length, interval);

          #ifdef DEBUG
          //std::cout << "length = " << length << std::endl;

          if ( neighbor_it == matrix.end() ) {
            std::cout << "Did not find entry in matrix!" << std::endl;
          }
          if ( i < length ) {
            //std::cout << "Edge is too long" << std::endl;
          }
          #endif // DEBUG

          if ( (neighbor_it != matrix.end()) && (i >= length) ) {

            #ifdef FSDEBUG
            std::cout << "edge cost = " << current_neighbor.cost << std::endl;
            std::cout << "accumulated cost = " << neighbor_it->second[i - length].cost << std::endl;
            #endif // DEBUG

            float nowcost = neighbor_it->second[i - length].cost + current_neighbor.cost;
            if ( nowcost < bestcost ) {
              bestcost = nowcost;
              best_predecessor = current_neighbor.st;
              best_steplength = length;
            }
          }
        }
      }

      fptas_backtrace fb = {bestcost, best_predecessor, best_steplength};
      it->second.push_back(fb);

      #ifdef FSDEBUG
      std::cout << "Best cost to state";
      it->first->print();
      std::cout << " in length " << i << " is " << bestcost << std::endl;
      #endif // DEBUG

      //if solution find break;
      if ( ( currentState->goalCheck(goalState, std::numeric_limits<float>::max() ) ) && ( bestcost <= c ) ) {

        std::vector<SearchState*> opt_path;

        int nowlength = i;
        std::cout << "Find sol finds solution at i = " << i << std::endl;

        while(currentState != startState){
          opt_path.push_back(currentState);
          fptas_backtrace fb = matrix[currentState][nowlength];
          currentState = fb.predecessor;
          nowlength -= fb.steplength;
        }
        opt_path.push_back(startState);


        solution sol = {opt_path, upscale(i, interval), bestcost};
        std::vector<solution> solutions;
        solutions.push_back(sol);
        return solutions;
      }
    }
  }

  std::cout << "Find sol fails.\n";
  throw 0;
}

float SearchSpace::decide_interval(float lb, float eps, float path_size){

  return (lb * eps / path_size);
  //return 1.0;
}

void SearchSpace::clear_map()
{
  for( auto it = vir_state_map.begin(); it != vir_state_map.end(); it++){
    delete it->second;
  }
  for( auto it = state_map.begin(); it != state_map.end(); it++){
    delete it->second;
  }
  state_map.clear();

  vir_state_map.clear();
  //no need to delete virtual states since they are deleted in state_map
}

void SearchSpace::reset_status(bool virtual_flag)
{
  if (virtual_flag) {
    for( auto it = vir_state_map.begin(); it != vir_state_map.end(); it++){
      it->second->reset_status();
    }    
  }
  else {
    for( auto it = state_map.begin(); it != state_map.end(); it++){
      it->second->reset_status();
    }     
  }
}

void SearchSpace::reset_values(bool virtual_flag)
{
  if (virtual_flag) {
    for( auto it = vir_state_map.begin(); it != vir_state_map.end(); it++){
      it->second->reset_values();
    }    
  }
  else {
    for( auto it = state_map.begin(); it != state_map.end(); it++){
      it->second->reset_values();
    }     
  }
}

solution SearchSpace::Larac_Advised_Astar(SearchState* startState, SearchState* goalState, float c, float eps, bool filter, bool user_study)
{
  clock_t t = clock();
  float elapsed_time = 0;

  reset_status(); // clean up all pruned status
  reset_values(); // clear the values

  //std::cout << "Start Running Larac Advised Astar...\n";
  solution astar_sol;
  astar_sol.path_cost = std::numeric_limits<float>::max();

  if (startState == NULL){ //check the inputs
    std::cout << "Error: Start State not in State Map\n";
    return astar_sol;
  }
  else if ( goalState == NULL ){ //check the inputs
    std::cout << "Error: Goal State not in State Map\n";
    return astar_sol;
  }

  elapsed_time = float(clock() - t)/float(CLOCKS_PER_SEC);
  if (elapsed_time > 10) { throw 1; }

  // try beta larac first to find whether greedy is enough
  larac_result l = larac_search(getLeader(startState), getLeader(goalState), c, eps, false);
  #ifdef DEBUG
  std::cout << "---LARAC results (Beta)---\n";
  std::cout << "Lambda: " << l.lambda << " , UB: " << l.sol.path_length << " , LB: " << l.obj_lower_bound << "\n";
  #endif // DEBUG

  bool beta_inf = (l.lambda == -2);
  bool beta_greedy = (l.lambda == -1);

  float beta_lambda = l.lambda; // remember lambda

  elapsed_time = float(clock() - t)/float(CLOCKS_PER_SEC);
  if (elapsed_time > 10) { throw 1; }

  if (l.lambda == -1) {
    std::cout << "This is a Beta-Greedy Instance.\n";
    l.lambda = 0;
    beta_greedy = true;
    astar_sol = l.sol;
  }

  else {
    l = larac_search(getLeader(startState), getLeader(goalState), c, eps, true); // try alpha to find if no solution exists
    #ifdef DEBUG
    std::cout << "---LARAC results (Alpha)---\n";
    std::cout << "Lambda: " << l.lambda << " , UB: " << l.sol.path_length << " , LB: " << l.obj_lower_bound << "\n";
    #endif

    if (l.lambda == -2) {
      std::cout << "This instance is infeasible!\n";
      if (filter) { throw 1; }
      return astar_sol;
    }

    else {
      if (l.lambda == -1) {
        //#ifdef DEBUG
        std::cout << "This is a Alpha-Greedy Instance.\n";
        //#endif
        l.lambda = 0;
      }

      astar_sol = heuristic_search(
      startState,
      goalState,
      true, //heuristic_flag
      false, //virtual_flag
      false, //cost_flag
      l.lambda, // A* using lambda from alpha LARAC.
      c,
      true // alpha_or_beta is not relevant here
      );

      elapsed_time = float(clock() - t)/float(CLOCKS_PER_SEC);
      if (elapsed_time > 10) { throw 1; }

      if (astar_sol.path_cost > c){ // alpha is not feasible
        std::cout << "Alpha Solution is not feasible.\n";
        reset_status(); // clean up all pruned status
        reset_values(); // clear the values
        reset_status(true);
        if (beta_greedy) {beta_lambda = 0;}
        astar_sol = heuristic_search(
        startState,
        goalState,
        true, //heuristic_flag
        false, //virtual_flag
        beta_inf, //cost_flag
        beta_lambda, // A* using lambda from beta LARAC.
        c,
        true // alpha_or_beta is not relevant here
        );
        elapsed_time = float(clock() - t)/float(CLOCKS_PER_SEC);
        if (elapsed_time > 10) { throw 1; }
      }

      if (astar_sol.path_cost > c){ // beta is still not feasible
        std::cout << "Beta Solution is not feasible. Surrender\n";
      }
    }
  }

  std::cout << "A* current solution: " << astar_sol.path_length << ", cost = " << astar_sol.path_cost << std::endl;

  if (user_study){
    PhyPath path = optimize_solution(astar_sol, startState); //Try MCPP
    astar_sol.path_cost = path.cost; // update the solution
    std::cout << "Optimized A* cost = " << astar_sol.path_cost << std::endl;
  }

  #ifdef DEBUG
  for (auto it = astar_sol.path.begin(); it != astar_sol.path.end(); ++it){
    (*it)->print();
  }
  #endif //DEBUG

  return astar_sol;
}

solution SearchSpace::cost_first(SearchState* startState, SearchState* goalState, float c, float eps, bool pruning_flag, bool lambda_flag, bool user_study)
{
  reset_status(); // clean up all pruned status
  reset_values(); // clear the values
  reset_status(true);

  solution cost_sol;
  cost_sol.path_cost = std::numeric_limits<float>::max();

  if (startState == NULL){ //check the inputs
    std::cout << "Error: Start State not in State Map\n";
    return cost_sol;
  }
  else if ( goalState == NULL ){ //check the inputs
    std::cout << "Error: Goal State not in State Map\n";
    return cost_sol;
  }

  cost_sol = heuristic_search(
    startState,
    goalState,
    false, //heuristic_flag
    false, //virtual_flag
    true, //cost_flag
    0, // lambda = 0
    0, // no cost constraint
    true // alpha_or_beta
    );

  std::cout << "CostFirst cost = " << cost_sol.path_cost << ", length = " << cost_sol.path_length << std::endl;

  if (user_study){
    PhyPath path = optimize_solution(cost_sol, startState); //Try MCPP
    cost_sol.path_cost = path.cost; // update the solution
    std::cout << "Optimized CostFirst cost = " << cost_sol.path_cost << std::endl;
  }

  return cost_sol;

}

solution SearchSpace::pruning_and_fptas(SearchState* startState, SearchState* goalState, float c, float eps, bool pruning_flag, bool lambda_flag, bool user_study)
{
  reset_status(); // clean up all pruned status
  reset_values(); // clear the values

  solution final_sol;
  solution cost_sol;
  final_sol.path_cost = std::numeric_limits<float>::max();
  int findsol_lb = 0;

  if (startState == NULL){ //check the inputs
    std::cout << "Error: Start State not in State Map\n";
    return final_sol;
  }
  else if ( goalState == NULL ){ //check the inputs
    std::cout << "Error: Goal State not in State Map\n";
    return final_sol;
  }

  float lambda, lb, ub, interval;

  // try beta larac first to find whether greedy is enough
  larac_result l = larac_search(getLeader(startState), getLeader(goalState), c, eps, false);

  #ifdef DEBUGG
  std::cout << "---LARAC results (Beta)---\n";
  std::cout << "Lambda: " << l.lambda << " , UB: " << l.sol.path_length << " , LB: " << l.obj_lower_bound << "\n";
  #endif // DEBUG

  if (l.lambda == -1) {
    std::cout << "This is a Beta-Greedy Instance.\n";
    l.lambda = 0;
    if (user_study){
      PhyPath path = optimize_solution(l.sol, startState); //MCPP on Larac Alpha
      l.sol.path_cost = path.cost; // update the solution
      std::cout << "Optimized Greedy cost = " << l.sol.path_cost << std::endl;
    }
    if (pruning_flag){
      return l.sol;
    }
  }

  //if(false){;}

  //else {
  // try CostFirst to find if no solution exists
    reset_status(); // clean up all pruned status
    reset_values(); // clear the values
    reset_status(true);
    cost_sol = heuristic_search(
    startState,
    goalState,
    false, //heuristic_flag
    false, //virtual_flag
    true, //cost_flag
    0, // lambda = 0
    0, // no cost constraint
    true // alpha_or_beta
    );

    findsol_lb = cost_sol.path_length; // this is the lower bound of obj, found by CostFirst
    interval = decide_interval(findsol_lb, eps, cost_sol.path.size()); // DP interval decided
    //if (!pruning_flag) {interval = 1;}
    findsol_lb = findsol_lb/interval + cost_sol.path.size(); // Add the error margins for max iteration of DP loops

    reset_status(); // clean up all pruned status
    reset_values(); // clear the values
    reset_status(true);

    l = larac_search(getLeader(startState), getLeader(goalState), c, eps, true); 

    #ifdef DEBUGG
    std::cout << "---LARAC results (Alpha)---\n";
    std::cout << "Lambda: " << l.lambda << " , UB: " << l.sol.path_length << " , LB: " << l.obj_lower_bound << "\n";
    #endif // DEBUGG

    if (cost_sol.path_cost > c) {
      std::cout << "This instance is infeasible!\n";
      return final_sol;
    }

    else {
      if (l.lambda == -1) {
        #ifdef DEBUGG
        std::cout << "This is an Alpha-Greedy Instance.\n";
        #endif
        l.lambda = 0;
      }

      lambda = l.lambda;
      lb = l.obj_lower_bound;
      ub = l.sol.path_length;

      #ifdef DEBUGG
      std::cout << "(Lambda) Heuristic solution = " << ub << std::endl;
      #endif // DEBUGG

      if (pruning_flag){
        ub = prune_space(startState, goalState, c, lambda_flag, lambda);
        #ifdef DEBUGG
        std::cout << "Pruning solution = " << ub << std::endl;
        #endif // DEBUGG
      }

      #ifdef DEBUGG
      std::cout << "Interval = " << interval << std::endl;
      #endif // DEBUGG
    }
  //}

  try{
    final_sol = find_sol(startState, goalState, c, interval, !pruning_flag, findsol_lb)[0];

    std::cout << "FPTAS current solution: " << final_sol.path_length << ", cost = " << final_sol.path_cost << std::endl;

    if (user_study){
      PhyPath path = optimize_solution(final_sol, startState);
      final_sol.path_cost = path.cost;
      std::cout << "Optimized cost = " << path.cost << std::endl;
    }

    //calculate the real path distance
    float real_dist = 0;
    int last_vx = 0;
    int last_vy = 0;
    for (auto it = final_sol.path.begin(); it != final_sol.path.end(); ++it){
      if (it != final_sol.path.begin()){
        float dx = (*it)->getGammaVX() - last_vx;
        float dy = (*it)->getGammaVY() - last_vy;
        real_dist += std::pow( (std::pow(dx,2.0)+std::pow(dy,2.0) ), 0.5 );
        std::cout << "adding" << std::pow( (std::pow(dx,2.0)+std::pow(dy,2.0) ), 0.5 ) << "to real_dist\n";
      }
      last_vx = (*it)->getGammaVX();
      last_vy = (*it)->getGammaVY();
      #ifdef DEBUGG
      (*it)->print();
      #endif // DEBUG
    }
    final_sol.path_length = real_dist;  
    return final_sol;
  }
  catch(int i){
    if(i==99){
      throw 99;
    }

    //if find_sol does not return a solution, return cost first
    return cost_sol;
  }
}

SearchState* SearchSpace::getState(loco_info lc)
{
  auto it = state_map.find(lc);
  if ( !(it == state_map.end() ) ) { return it->second; }
  else {
    std::cout << "Error: state not found.\n";
    return NULL;
  }
}

SearchState* SearchSpace::getVirState(loco_info lc)
{
  auto it = vir_state_map.find(lc);
  if ( !(it == vir_state_map.end() ) ) { return it->second; }
  else {
    std::cout << "Error: state not found.\n";
    return NULL;
  }
}

PhyPath SearchSpace::optimize_solution(solution& sol, SearchState* startState)
{
  PhyPath dummy_path;
  if (sol.path.size() == 0) {
    std::cout << "Error: Cannot optimize an empty solution.\n";
    return PhyPath();
  }
  std::reverse(sol.path.begin(), sol.path.end()); // reverse the solution path

  for (unsigned i = 0; i < sol.path.size(); i++){
    std::cout << "(" << sol.path[i]->getGammaVX() << "," << sol.path[i]->getGammaVY() << ")->";
  }

  std::cout << std::endl;
  float last_angle = startState->getThetaV();
  std::map<Point, float> vir_turning_angles;

  int dx, dy;
  std::vector<std::pair<float, float> > vsteps;
  for (unsigned i = 0; i < sol.path.size() - 1; i++){
    dx = sol.path[i+1]->getGammaVX() - sol.path[i]->getGammaVX();
    dy = sol.path[i+1]->getGammaVY() - sol.path[i]->getGammaVY();

    float l = std::pow( std::pow(dx, 2.0) + std::pow(dy, 2.0), 0.5);
    float temp_angle = (std::atan2(dy, dx)/M_PI)*180.;
    float angle = (temp_angle<0)?temp_angle+360.:temp_angle;

    float angle_diff = angle - last_angle;
    if (angle_diff < -180.) {angle_diff += 360.;}
    else if (angle_diff > 180.) {angle_diff -= 360.;}
    Point p(sol.path[i]->getGammaVX(), sol.path[i]->getGammaVY());
    vir_turning_angles[p] = angle_diff;
    last_angle = angle;

    vsteps.push_back(std::make_pair(l, angle));
  }

  Point startP(startState->getGammaPX(), startState->getGammaPY());
  float phy_theta = startState->getThetaP();
  float vir_theta = startState->getThetaV();

  PhyPath path;
  float final_cost = 0;

  #ifdef DEBUG_MCPP
  std::cout << "V Steps:\n";
  for (auto vv = vsteps.begin(); vv != vsteps.end(); vv++){
    std::cout << "l = " << vv->first << ", theta = " << vv->second << "\n";
  }
  std::cout << "start P = (" << startP.x << "," << startP.y << ")\n phy_theta = " << phy_theta << "\n vir_theta = " << vir_theta
  << "\n";
  #endif
  //std::cout << "call greedy_optimize_path" << std::endl;

  Point startV(startState->getGammaVX(), startState->getGammaVY());
  final_cost = p_ptr->Greedy_Optimize_Path(startP, phy_theta, vsteps, vir_theta, startV, vir_turning_angles);

  dummy_path.cost = final_cost;
  return dummy_path;
}
