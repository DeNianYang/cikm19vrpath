#include "VirIndex.h"
#include "PhyIndex.h"
using namespace std;

int main() {

    LocoState init_state(6, 14, 0, 0, 0, 0);
    Point goal = Point(10, 2);
    DPCvsp cvsp(init_state, goal, 100);
    cvsp.read_phy_map("./map/phy.map");
    cvsp.read_vir_map("./map/supermarket.map");

#ifdef DEBUG
    cvsp.print_map("phy");
    cvsp.print_map("vir");
#endif // DEBUG

    cvsp.DP_CVSP();
    cvsp.print_path();
    cout << "\n";
    VirIndex vir_ind = VirIndex();

    return 0;
}
