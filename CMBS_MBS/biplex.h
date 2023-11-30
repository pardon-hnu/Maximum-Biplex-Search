#include<iostream>
#include<list>
#include<vector>
#include<time.h>
using namespace std;

class biplexminer{

    public:
        int res_count;
        int res_size;
        void Miner();
        biplexminer(int **Graph, int Graph_size, int Bipartite_index,int *degree, int S, int k);
        ~biplexminer();
        vector<int> Lnode;
        vector<int> Rnode;
        vector<vector<int>> Xs;
        double GenerateT;
        double EnumerateT;

    private:
        int **Graph;
        int Graph_size;
        int Bipartite_index;
        int S;
        int k;
        int *G_State;
        int *G_index;
        int *G_label;
        int *G_temp;
        int *G_prune;
        int *degree;

        void Enumerate(vector<int> &X, vector<int> &N_X, vector<int> &Cand_p, vector<int> &Cand_q);
        void Generate(vector<int> &X, vector<int> &N_X,vector<int> &Cand_p, vector<int> &Cand_q);
        void List(vector<int> &X, vector<int> var_Y, vector<int> &con_Y, vector<int> &cur_Y,vector<int> &candY_q, vector<int> &ext_p, vector<int> &ext_q);
        bool CheckMaximality(vector<int> &X, vector<int> &con_Y, vector<int> &cur_Y, vector<int> &candY_q, vector<int> &ext_p, vector<int> &ext_q);
        void Special(vector<int> &X, vector<int> &var_Y, vector<int> &con_Y, vector<int> &cur_Y,vector<int> &candY_q, vector<int> &ext_p, vector<int> &ext_q);
	    void Output(vector<int> &L, vector<int> &R1, vector<int> &R2);
	    int GetSize(vector<int> &L, vector<int> &R1, vector<int> &R2);
};
