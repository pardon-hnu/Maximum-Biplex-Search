#include<iostream>
#include<list>
#include<map>
#include<algorithm>
#include <vector>
using namespace std;
struct xypair
{
    private:
        int x;
        int y;
        int ub;
        int lb;
    public:
        xypair()
        {
            this->x=0;
            this->y=0;
            this->ub=0;
            this->lb=0;
        }
        xypair(int x,int y,int k)
        {
            this->x=x;
            this->y=y;
            this->ub=x*y;
            this->lb=(x>y)?(x*y-y*k):(x*y-x*k);
        }
        int getx()
        {
            return this->x;
        } 
        int gety()
        {
            return this->y;
        }
        int getupperbound()
        {
            return this->ub;
        }
        int getlowerbound()  
        {
            return this->lb;
        }
};
class Core{
    public:
        Core(int **Graph, int *degree, int graph_size, int bipartite, int K);
        int CorePrune(int **&new_graph, int *&degree, int &new_bi, int core_number);
        int XYCorePrune(int **&new_graph, int *&degree, int &new_bi, int alphi, int beta);
        void Coredecompose();
        int GetMaxcore();
        vector<xypair> GetPairs(int left,int right);
        int GetMaxY(int x);
        int GetMaxX(int y);
        int *G_index;
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
Core::Core(int **Graph, int *degree, int graph_size, int bipartite, int K){
    this->Graph=Graph;
    this->degree=degree;
    this->graph_size=graph_size;
    this->K=K;
    this->bipartite=bipartite;

    int max=0;
    int Umax = 0;
    int Vmax = 0;
    G_temp=new int[graph_size];
    G_label=new int[graph_size];
    G_index=new int[graph_size];
    G_state=new int[graph_size];

    for(int i=0;i<bipartite;++i){
        G_state[i]=1;
        G_index[i]=-1;
        G_label[i]=0;
        G_temp[i]=degree[i];
        if(degree[i]>Umax){
            Umax=degree[i];
        }
        if(degree[i]>max){
            max=degree[i];
        }
    }
    for(int i=bipartite;i<graph_size;++i){
        G_state[i]=1;
        G_index[i]=-1;
        G_label[i]=0;
        G_temp[i]=degree[i];
        if(degree[i]>Vmax){
            Vmax=degree[i];
        }
        if(degree[i]>max){
            max=degree[i];
        }
    }
    this->max_u=Umax;
    this->max_v=Vmax;
    this->max_degree=max;
}

void Core::Coredecompose(){

    int *bin=new int[max_degree+1];
    for(int i=0;i<=max_degree;++i){
        bin[i]=0;
    }
    for(int i=0;i<graph_size;++i){
        bin[degree[i]]+=1;
    }


    int start=0;
    for(int d=0;d<=max_degree;++d){
        int num=bin[d];
        bin[d]=start;
        start+=num;
    }

    int *pos=new int[graph_size+1];
    int *vert=new int[graph_size+1];
    for(int i=0;i<graph_size+1;++i){
        pos[i]=0;
        vert[i]=0;
    }
    for(int i=0;i<graph_size;++i){
        pos[i]=bin[degree[i]];
        vert[pos[i]]=i;
        bin[degree[i]]+=1;
    }

    for(int i=max_degree;i>=1;--i){
        bin[i]=bin[i-1];
    }
    bin[0]=1;


    for(int i=0;i<graph_size;++i){
        int node=vert[i];
        G_label[node]=G_temp[node];
        for(int j=0;j<degree[node];++j){
            int neg=Graph[node][j];
            if(G_temp[neg]>G_temp[node]){
                int du=G_temp[neg]; int pu=pos[neg];
                int pw=bin[du]; int w=vert[pw];
                if(neg!=w){
                    pos[neg]=pw; vert[pu]=w;
                    pos[w]=pu; vert[pw]=neg;
                }
                bin[du]+=1;
                G_temp[neg]-=1;

            }
        }
    }

}

int Core::GetMaxcore(){
//    for(int i=0;i<graph_size;i++){
//        printf("node: %d, coreness: %d, degree: %d.\n", i, G_temp[i], degree[i]);
//    }
    int max=0;
    for(int i=0;i<graph_size;++i){
        if(G_label[i]>max){
            max=G_label[i];
        }
    }

    int temp_count=0;
    for(int i=0;i<graph_size;++i){
        if(G_label[i]>=max){
            temp_count++;
        }
    }

    // cout<<1.0*temp_count/graph_size<<endl;
    // cout<<graph_size<<endl;
    // cout<<temp_count<<endl;

    return max;
}

int Core::CorePrune(int **&new_graph, int *&new_degree, int &new_bi, int core_number){

    int count=0;
    int bi_count=0;

    for(int i=0;i<graph_size;++i){
        if(G_label[i]>=core_number){
            G_index[i]=count;
            count++;
            if(i<bipartite){
                bi_count++;
            }
        }else{
            G_index[i]=-1;
        }
    }

    new_bi=bi_count;
    new_degree=new int[count];


    new_graph=new int*[count];
    for(int i=0;i<graph_size;++i){
        if(G_index[i]>=0){
            int temp_count=0;
            for(int j=0;j<degree[i];++j){
                if(G_index[Graph[i][j]]>=0){
                    temp_count++;
                }
            }
            int *neg=new int[temp_count];
            new_degree[G_index[i]]=temp_count;
            temp_count=0;

            for(int j=0;j<degree[i];++j){
                if(G_index[Graph[i][j]]>=0){
                    neg[temp_count]=G_index[Graph[i][j]];
                    temp_count++;
                }
            }
            new_graph[G_index[i]]=neg;
        }
    }

    return count;
}


int Core::XYCorePrune(int **&new_graph, int *&new_degree, int &new_bi, int alphi, int beta){

    for(int i=0;i<graph_size;++i){
        G_temp[i]=degree[i];
        G_state[i] = 1;
    }

    //xycore
    while (1){
        int flag1 = 1;
        int flag2 = 1;
        for (int i = 0; i < bipartite; i++){
            if (G_state[i]>0){
                if(G_temp[i]<alphi){
                    flag1 = 0;
                    G_state[i] = 0;
                    for (int j = 0; j < degree[i]; j++){
                        G_temp[Graph[i][j]] -= 1;
                    }
                }
            }
        }
        for (int i = bipartite; i < graph_size; i++){
            if (G_state[i]>0){
                if(G_temp[i]<beta){
                    flag2 = 0;
                    G_state[i] = 0;
                    for (int j = 0; j < degree[i]; j++){
                        G_temp[Graph[i][j]] -= 1;
                    }
                }
            }
        }
        if (flag1&&flag2) break;
    }

    int count = 0;
    int bi_count = 0;

    for (int i=0;i<graph_size;i++){
        if (G_state[i] >0){
            G_index[i] = count;
            count +=1;
            if(i<bipartite) {bi_count++;}
        }else{
            G_index[i] = -1;
        }
        //cout<<"i: "<<i<<" index: "<<G_index[i]<<endl;
    }

    new_graph = new int*[count];
    new_bi = bi_count;
    new_degree = new int[count];

    for(int i=0;i<graph_size;i++){
        if (G_index[i] >=0){
            int temp_count = 0;
            for (int j=0;j<degree[i];j++){
                if (G_index[Graph[i][j]]>=0)
                    temp_count++;
            }
            //���ھӷ���������
            int *neg = new int[temp_count];
            new_degree[G_index[i]] = temp_count;

            //���ھ���������µ�graph��
            temp_count = 0;
            for (int j=0; j<degree[i]; j++){
                if (G_index[Graph[i][j]] >=0){
                    neg[temp_count] = G_index[Graph[i][j]];
                    temp_count++;
                }
            }
            new_graph[G_index[i]] = neg;
        }
    }
    //cout<<"graphsize: "<<count<<" bipartite: "<<new_bi<<endl;
    int size = 0;
    for(int i=0;i<count;i++)
        size+=new_degree[i];
    // cout<<"U: "<<bi_count<<", V: "<<count-bi_count<<", core size:"<<size/2<<endl;
    return count;
}

vector<xypair> Core::GetPairs(int left,int right){
    vector<xypair> pairs;
    for(int x = left; x<max_v+K; x++){
        int max_y = GetMaxY(x);
        if(max_y < 2*K+1)
            break;
        // cout<<"now x: "<<x<<" and max y: "<<max_y<<endl;
        for(int y=right;y<=max_y;y++)
        {
            xypair p(x,y,K);
            pairs.push_back(p);
        }
        
    }
    return pairs;
}


int Core::GetMaxY(int x){
    int temp_y = 2*K+1;
    int num_u = bipartite;
    int flag1 = 1;
    int flag2 = 1;

    for(int i=0;i<graph_size;++i){
        G_temp[i]=degree[i];
        G_state[i] = 1;
    }

    while (num_u){
        flag1 = 1;
        flag2 = 1;

        for (int i = 0; i < bipartite; i++){
            if (G_state[i]>0){
                if(G_temp[i]<temp_y-K){
                    flag1 = 0;
                    num_u --;
                    G_state[i] = 0;
                    for (int j = 0; j < degree[i]; j++){
                        G_temp[Graph[i][j]] -= 1;
                    }
                }
            }
        }

        for (int i = bipartite; i < graph_size; i++){
            if (G_state[i]>0){
                if(G_temp[i]<x-K){
                    flag2 = 0;
                    G_state[i] = 0;
                    for (int j = 0; j < degree[i]; j++){
                        G_temp[Graph[i][j]] -= 1;
                    }
                }
            }
        }
        if (flag1&&flag2) temp_y++;
    }
    //pairs[x] = temp_y;
    //cout<<"x: "<<x<<"y: "<<temp_y<<endl;
    return temp_y-1;
}

int Core::GetMaxX(int y){
    int temp_x = 2*K+1;
    int num_u = graph_size-bipartite;
    int flag1 = 1;
    int flag2 = 1;

    //重置degree和节点状态
    for(int i=0;i<graph_size;++i){
        G_temp[i]=degree[i];
        G_state[i] = 1;
    }



    while (num_u){
        flag1 = 1;
        flag2 = 1;
        //移除Y中的节点
        for (int i = bipartite; i < graph_size; i++){
            if (G_state[i]>0){
                if(G_temp[i]<temp_x-K){
                    flag1 = 0;
                    num_u --;
                    G_state[i] = 0;
                    for (int j = 0; j < degree[i]; j++){
                        G_temp[Graph[i][j]] -= 1;
                    }
                }
            }
        }

        //移除X中的节点
        for (int i = 0; i < bipartite; i++){
            if (G_state[i]>0){
                if(G_temp[i]<y-K){
                    flag2 = 0;
                    G_state[i] = 0;
                    for (int j = 0; j < degree[i]; j++){
                        G_temp[Graph[i][j]] -= 1;
                    }
                }
            }
        }
        if (flag1&&flag2) temp_x++;
    }
    //pairs[x] = temp_y;
    // cout<<"x: "<<temp_x<<"y: "<<y<<endl;
    return temp_x-1;
}
