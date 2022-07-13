#include "PhyIndex.h"
using namespace std;

int main() {
    cout << "Start Running\n";
    PhyIndex phy_ind = PhyIndex(5);
    phy_ind.read_graph("./map/phy.map");
    phy_ind.print_graph();
    phy_ind.build_phy_index();
    phy_ind.export_to_file("phy.index");

    /*DPCvsp cvsp = DPCvsp();
    cvsp.read_phy_map("./map/phy.map");
    cvsp.print_map("phy");*/
    return 0;
}
