#include<iostream>
#include<list>
#include<vector>

using namespace std;

class biplexminer_xy{

    public:
        int res_count;
        int res_size;
        void Miner();
        void Miner1();
        biplexminer_xy(int **Graph, int Graph_size, int Bipartite_index,int *degree, int alpha, int beta, int k);
        ~biplexminer_xy();
        vector<int> maxL;
        vector<int> maxR;
    private:
        int **Graph;
        int Graph_size;
        int Bipartite_index;
        int alpha;
        int beta;
        int k;
        int *G_State;
        int *G_index;
        int *G_label;
        int *G_temp;
        int *G_prune;
        int *degree;


        void Enumerate(vector<int> &X, vector<int> &N_X, vector<int> &Cand_p, vector<int> &Cand_q);
        void Enumerate1(vector<int> &X, vector<int> &N_X, vector<int> &Cand_p, vector<int> &Cand_q);
        void Generate(vector<int> &X, vector<int> &N_X,vector<int> &Cand_p, vector<int> &Cand_q);
        void List(vector<int> &X, vector<int> var_Y, vector<int> &con_Y, vector<int> &cur_Y,vector<int> &candY_q, vector<int> &ext_p, vector<int> &ext_q);
        bool CheckMaximality(vector<int> &X, vector<int> &con_Y, vector<int> &cur_Y, vector<int> &candY_q, vector<int> &ext_p, vector<int> &ext_q);
        void Special(vector<int> &X, vector<int> &var_Y, vector<int> &con_Y, vector<int> &cur_Y,vector<int> &candY_q, vector<int> &ext_p, vector<int> &ext_q);
        int GetSize(vector<int> &L, vector<int> &R1, vector<int> &R2);
        void Output(vector<int> &L, vector<int> &R1, vector<int> &R2);
};

