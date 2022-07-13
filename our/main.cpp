#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <string>
#include <thread>
#include "SearchState.h"
#include "LocoState.h"
#include "SearchSpace.h"
#include "PhyIndex.h"
#include "kshortestpath.h"
#include "VirIndex.h"
#include "Tools.h"

using namespace std;

unsigned max_lv = 0;
unsigned samples = 0;
unsigned clb = 0;
unsigned cub = 0;
unsigned k = 5;
bool build_vir = false;
bool build_phy = false;
bool simplify = true; //for now, all implements are simplified
bool query = true;
bool filter = false;
bool user_study = false;
int query_lb = 1;
int query_ub = 10000;
string outfilename = "output";
bool naiveDP = false;
int KSPflag = 0;

const int LANA = 0;
const int LANA_ADV = 1;
const int NAIVE_DP = 2;
const int COST = 3;
const int K_SP = 4;

string method_names [5] = {"LANA", "LANA-ADV", "NAIVE_DP", "COST", "KSP"};
float compute_time [5] = {};
float objective [5] = {};
float feasible_objective [5] = {};
float rw_cost [5] = {};
int feasible_times [5] = {};
int return_sol_times [5] = {};

int ALL_feasible = 0;
int Incomplete_Queries = 0;
int Valid_Queries = 0;

float RESET_COST = 5.0;

void print_results(ofstream& of){
	of << endl;

	of << "=============================================================================================================" << endl;

	of << "Run Naive DP: " << ((naiveDP) ? "Y" : "N ");
	of << ", Run KSP: " << ((KSPflag > 0) ? "Y" : "N ");
	of << ", Filter Impossible: " << ((filter > 0) ? "Y" : "N ");
	of << ", Reset Cost Ratio: " << RESET_COST;

	/*if (naiveDP){ 		of << "Y";
	}
	else{
		of << "N";
	}

	of << ", Run KSP: ";

	if (KSPflag > 0){
		of << "Y";
	}
	else{
		of << "N";
	}*/

	of << endl;
	of << "Successful Queries = " << Valid_Queries << ", All Feasible Queries = " << ALL_feasible << endl;

	of << "==METHOD== | == AVG. == | == FEAS == | == FEAS == | == AVG. == | ==ALLFEA== | == AVG. == " << endl;
	of << "== NAME == | == TIME == | == TIMES== | == RATIO== | == OBJ. == | == OBJ. == | == COST == " << endl;
	of << "-----------|------------|------------|------------|------------|------------|------------" << endl;
	for(int i = 0; i < 5; i++){
		if (((i==2) && (!naiveDP)) || ((i==4) && (KSPflag == 0))){continue;}
		of << setw(10) << method_names[i];
		of << " | " << setprecision(4) << setw(10) << compute_time[i]/float(feasible_times[3]);
		of << " | " << setprecision(0) << setw(10) << feasible_times[i];
		of << " | " << setprecision(2) << setw(9)  << 100 * float(feasible_times[i])/float(feasible_times[3]) << '%';
		of << " | " << setprecision(3) << setw(10) << objective[i]/float(feasible_times[i]);
		of << " | " << setprecision(3) << setw(10) << feasible_objective[i]/float(ALL_feasible);
		of << " | " << setprecision(3) << setw(10) << rw_cost[i]/float(return_sol_times[i]);
		of << endl;
	}
	of << "=============================================================================================================" << endl;
}

int main(int argc, char const *argv[])
{
	if (argc < 3) {
		cout << "Error: missing map files.\n";
		return -1;
	}

	const char* virtualmap = argv[1];
	const char* physicalmap = argv[2];

	for (int i = 2; i < argc; i++){
		const char* arg = argv[i];
		if (arg[0] == '-'){
			switch (arg[1]){
				case 'n':{
					i++;
					samples = atoi(argv[i]);
					break;
				}
				case 'l':{
					i++;
					max_lv = atoi(argv[i]);
					break;
				}
				case 'c':{
					i++;
					clb = atoi(argv[i]);
					i++;
					cub = atoi(argv[i]);
					break;
				}
				case 'd':{
					naiveDP = true;
					break;
				}
				case 'o':{
					i++;
					outfilename = argv[i];
					break;
				}
				case 's':{
					simplify = true;
					break;
				}
				case 'u':{
					user_study = true;
					break;
				}
				case 'f':{
					filter = true;
					break;
				}
				case 'i':{
					query = false;
					break;
				}
				case 'v':{
					build_vir = true;
					break;
				}
				case 'p':{
					build_phy = true;
					break;
				}
				case 'k':{
					i++;
					KSPflag = atoi(argv[i]);
					break;
				}
				case 'r':{
					i++;
					RESET_COST = atoi(argv[i]);
					break;
				}
				case 'q':{
					i++;
					query_lb = atoi(argv[i]);
					i++;
					query_ub = atoi(argv[i]);
					break;
				}
				default:
				;
			}
		}
	}

	string vindexname = outfilename;
	vindexname.append("_vindex");
	string pindexname = outfilename;
	pindexname.append("_pindex");
	string kspvinputname = outfilename;
	kspvinputname.append("_ksp_input");
	string kspindexname = outfilename;
	kspindexname.append("_ksp_index");
	string colafilename = outfilename;
	colafilename.append("_cola_queries");
	string userstudymapname = outfilename;
	userstudymapname.append("_user_study_map");
	string userstudypathname = outfilename;
	userstudypathname.append("_user_study_path");


	if (!query) {
		samples = 0;
	}
	else if (max_lv * samples * cub == 0) {
		cout << "Error: Missing numerical parameters.\n";
		return -1;
	}

	if (clb > cub) {
		cout << "Error: lower bound is higher than upper bound.\n";
		return -1;
	}

	cout << "Virtual map: " << virtualmap << endl;
	cout << "Physical map: " << physicalmap << endl;

	ofstream of,cola;
    of.open(outfilename);
    cola.open(colafilename);
    cola << samples << endl;

	clock_t t = clock();

	VirIndex v(max_lv);
	v.read_graph(virtualmap);
	if (user_study){
		v.output_user_study_map(userstudymapname);
	}

	float max_query_length = pow(pow(v.width, 2)+pow(v.length, 2), 0.5);
	if (max_query_length < query_lb){
		of << "Query group out of bound.\n";
		return 1;
	}

	float vir_index_time = 0;

	if (build_vir){
		v.build_vir_index();
		v.export_to_file(vindexname);
		v.export_for_ksp(kspvinputname);
		vir_index_time = float(clock() - t)/CLOCKS_PER_SEC;
	    of << "Virtual Index Finished and Exported in " << float(clock() - t)/CLOCKS_PER_SEC << " seconds." << endl;
	    t = clock();

	    if(KSPflag > 1){
			try{
				KShortestPath ksp;
			    ksp.read_graph(kspvinputname.c_str(), k);
			    ksp.read_type(virtualmap);
			    ksp.construct_index();
			    ksp.write_index(kspindexname.c_str());
			    ksp.add_indexing_time(vir_index_time);
			    ksp.print_summary();

			    of << "K-Shortest-Path Index Finished and Exported in " << float(clock() - t)/CLOCKS_PER_SEC << " seconds." << endl;
			    cout << "K-Shortest-Path Index Finished and Exported.\n";
			    t = clock();
			}
			catch(...){
				cout << "K-Shortest-Path building index failure. Return." << endl;
				return 0;
			}
		}
	}

	t = clock();

	PhyIndex p(max_lv);
	p.read_graph(physicalmap);
	p.pathfilename = userstudypathname;

	if (build_phy){
		p.build_sh_index(simplify);
		p.export_sh_index(pindexname);
	    of << "ShaoHeng's Physical Index Finished and Exported in " << float(clock() - t)/CLOCKS_PER_SEC << " seconds." << endl;
	}

	t = clock();

	KShortestPath ksp;

	if(KSPflag > 0){
	    double query_time;
	    ksp.read_index(kspindexname.c_str());
	    of << "K-Shortest-Path Index Read in " << float(clock() - t)/CLOCKS_PER_SEC << " seconds." << endl;
	    cout << "K-Shortest-Path Index Read.\n";
	}

    t = clock();

	SearchSpace s(max_lv, p);
	bool flag = s.build_states(virtualmap, physicalmap, vindexname, pindexname, simplify);

	if(flag){
		of << "Search environment built in " << float(clock() - t)/CLOCKS_PER_SEC << " seconds." << endl;

		if (!query) { return 0; } // no query will be performed

		of << "Query Log:\n";

		srand(time(0));

		while (Valid_Queries < samples){
			int vx1, vy1, vx2, vy2;
			while(true){
				while(true){
					vx1 = rand() % v.width;
					vy1 = rand() % v.length;
					if ( v.is_valid_input(vx1, vy1) ) { break; }
				}

				while(true){
					vx2 = rand() % v.width;
					vy2 = rand() % v.length;
					if ( (vx1 == vx2) && (vy1 == vy2)) { continue; }
					if ( v.is_valid_input(vx2, vy2) ) { break; }
				}
				float dx = abs(vx1 - vx2);
				float dy = abs(vy1 - vy2);
				float dd = pow(pow(dx, 2) + pow(dy,2), 0.5);
				if ( (dd < query_lb) || (dd > query_ub) ) { continue; }
				break;
			}
			int px, py;
			int ptheta = 0;
			if (user_study){ // hardcoding for user study
				px = 2;
				py = 2;
			}
			else{
				while(true){
					px = rand() % p.width;
					py = rand() % p.length;
					if (!simplify) { ptheta = (rand() % 8) * 45; } //angle = 0 if simplified verision
					if ( p.is_valid_input(px, py) ) { break; }
				}
			}

			loco_info lc_s = {vx1, vy1, 0, px, py, ptheta};
			loco_info lc_t = {vx2, vy2, 0, px, py, ptheta}; //only vx2 and vy2 matter

			int c = clb + ( rand() % (cub-clb) );

			cout << endl << "Sampling start state = (" << vx1 << ',' << vy1 << ",0, " << px << ',' << py << ',' << ptheta << ")\n";
			cout << "Sampling destination position = (" << vx2 << ',' << vy2 << "), c = " << c << "\n";

			float temp_time [5] = {};
			solution temp_sol [5];
			for (int i=0; i<5; i++){
				temp_sol[i].path_length = 999999;
				temp_sol[i].path_cost = 999999;
			}

			try{
				t = clock();
				temp_sol[LANA] = s.Larac_Advised_Astar(s.getState(lc_s), s.getState(lc_t), c, 0.1, filter, user_study);
				temp_time[LANA] = float(clock() - t)/CLOCKS_PER_SEC;
				cout << "LANA ended in " << temp_time[LANA] << " seconds." << endl << endl;
			}
			catch(int i){
				if (i==1){ // only if instance is definitely infeasible, and filter = true
					continue;
				}
			}
			catch(...){
				Incomplete_Queries++;
				continue;
			}

			try{
				t = clock();
				temp_sol[LANA_ADV] = s.pruning_and_fptas(s.getState(lc_s), s.getState(lc_t), c, 0.1, true, false, user_study);
				temp_time[LANA_ADV] = float(clock() - t)/CLOCKS_PER_SEC;
				cout << "LANA_ADV ended in " << temp_time[LANA_ADV] << " seconds." << endl << endl;
			}
			catch(...){
				Incomplete_Queries++;
				continue;
			}

			if (naiveDP){
				try{
					t = clock();
					temp_sol[NAIVE_DP] = s.pruning_and_fptas(s.getState(lc_s), s.getState(lc_t), c, 0.1, false, false, user_study); // no pruning before NaiveDP
					temp_time[NAIVE_DP] = float(clock() - t)/CLOCKS_PER_SEC;
					cout << "NaiveDP ended in " << temp_time[NAIVE_DP] <<  " seconds." << endl << endl;
				}
				catch(...){
					Incomplete_Queries++;
					continue;
				}
			}
			else{
				temp_sol[NAIVE_DP].path_cost = 0;
				temp_sol[NAIVE_DP].path_length = 0;
			}

			try{
				t = clock();
				temp_sol[COST] = s.cost_first(s.getState(lc_s), s.getState(lc_t), c, 0.1, false, false, user_study); // no pruning before NaiveDP
				temp_time[COST] = float(clock() - t)/CLOCKS_PER_SEC;
				cout << "Cost First ended in " << temp_time[COST] <<  " seconds." << endl << endl;
			}
			catch(...){
				Incomplete_Queries++;
				continue;
			}

			bool kspfound = false;
			float kspdistance = 0;
			float kspcost = 0;

			if (KSPflag > 0){
				try{
					t = clock();
					VirPathes k_pathes = ksp.query(vx1, vy1, vx2, vy2);
					float p_size = min(p.width, p.length);

				    Point s = {px, py};
				    for(int i = 0; i < k_pathes.size();i++){
				    	VirPath now_path = k_pathes[i];
				    	kspdistance = 0;
				    	kspcost = 0;
				    	for (int j = 0; j < now_path.size() - 1; j++){
				    		float dx = now_path[j+1].x - now_path[j].x;
				    		float dy = now_path[j+1].y - now_path[j].y;
				    		kspdistance += pow( pow(dx, 2.0) + pow(dy, 2.0), 0.5);
				    		kspcost += ceil(kspdistance / p_size) * RESET_COST;
				    	}
				    	if ( kspcost <= c ){
				    		kspfound = true;
				    		temp_sol[K_SP].path_cost = kspcost;
				    		temp_sol[K_SP].path_length = kspdistance;
				            break;
				        }
		   			}			
		   			temp_time[K_SP] = float(clock() - t)/float(CLOCKS_PER_SEC);
					cout << "K-Shortest-Path query ended in " << temp_time[K_SP] << " seconds." << endl << endl;
				}
				catch(...){
				Incomplete_Queries++;
				continue;
				}
			}
			else{
				temp_sol[K_SP].path_length = 0;
				temp_sol[K_SP].path_cost = 0;
			}


			//if no exceptions thrown, record
			Valid_Queries++;
			of << "(" << setw(3) << vx1 << ',' << setw(3) << vy1 << ")->(" << setw(3) << vx2 << ',' << setw(3) << vy2 << "), c = " << setw(3) << c << " | ";

			for(int i=0; i<5; i++){
				of << "T=" << fixed << setprecision(4) << temp_time[i];
				if (temp_sol[i].path_cost <= c) {
					feasible_times[i]++;
					of << " / fea / ";
					of << "L=" << temp_sol[i].path_length << ' ';
				}
				else{ of << " / inf / L=0.0000 "; }
				of << " | ";
			}
			of << endl;

			cola << vx1 << ' ' << vy1 << ' ' << vx2 << ' ' << vy2 << ' ' << c;

			bool all_fea = true; //flag

			for(int i=0; i<5; i++){
				if (temp_sol[3].path_cost <= c){ //if the instance is feasible, record computation time
					compute_time[i] += temp_time[i];
				}
				if (temp_sol[i].path_cost <= c){ //if the method is feasible, record feasible obj
					//objective[i] += temp_sol[i].path_length;
				}
				else { all_fea = false; }
				if (temp_sol[i].path_cost < 10000){ //not infinite
					return_sol_times[i]++;
					rw_cost[i] += temp_sol[i].path_cost;
				}
			}

			if (all_fea){
				ALL_feasible++;
				for(int i=0; i<5; i++){
					feasible_objective[i] += temp_sol[i].path_length;
				}
				cola << ' ' << ceil(temp_sol[0].path_length) << " 1" << endl;	
			}
			else{
				cola << " 0 0" << endl;
			}

			if (temp_sol[3].path_cost <= c){
				for(int i=0; i<5; i++){
					if (temp_sol[i].path_length > 10000) { continue; }
					objective[i] += (temp_sol[i].path_length)/(temp_sol[3].path_length);
				}
				cola << ' ' << ceil(temp_sol[0].path_length) << " 1" << endl;
			}

			if (Valid_Queries % 20 == 0){
				print_results(of);
			}
		}
	}
	else {
		of << "Search environment was not built sucessfully.\n";
		return 1;
	}
	print_results(of);
}



