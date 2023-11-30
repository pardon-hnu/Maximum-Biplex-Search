#include<iostream>
#include<list>
#include<map>
#include<algorithm>

using namespace std;

class Core{

    public:
        Core(int **Graph, int *degree, int graph_size, int bipartite, int K);
        int CorePrune(int **&new_graph, int *&degree, int &new_bi, int core_number);
        int XYCorePrune(int **&new_graph, int *&degree, int &new_bi, int alphi, int beta);
        void Coredecompose();
        int GetMaxcore();
        map<int,int> GetPairs();
        int GetMaxY(int x);
        int GetMaxX(int y);
        int *G_index;
        int *Re_index;
    private:
        int **Graph;
        int *degree;
        int *G_temp;
        int *G_label;
        int *G_state;
        int K;
        int graph_size;
        int max_degree;
        int bipartite;
        int max_u;
        int max_v;

};

