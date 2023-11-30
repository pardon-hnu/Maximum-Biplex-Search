#pragma once
#include<iostream>
#include<list>
#include<algorithm>

using namespace std;

class CoreLocate{

    public:
        CoreLocate(int **Graph, int *degree, int graph_size, int bipartite, int K);
        int CorePrune(int **&new_graph, int *&degree, int &new_bi);
        void Coredecompose();
        int GetMaxcore();
        void CoreOrdering(int **&new_graph, int *&new_degree);
        void Bipartite_CoreOrdering(int **&new_graph, int *&new_degree);
        void Bipartite_R_CoreOrdering(int **&new_graph, int *&new_degree);
        void GetReIndex();
        int *G_order;
        int *G_index;
        int *re_index;

    private:
        int **Graph;
        int *degree;
        int *G_temp;
        int *G_label;
        
        int K;
        int graph_size;
        int max_degree;
        int bipartite;

};

CoreLocate::CoreLocate(int **Graph, int *degree, int graph_size, int bipartite, int K){
    this->Graph=Graph;
    this->degree=degree;
    this->graph_size=graph_size;
    this->K=K;
    this->bipartite=bipartite;

    int max=0;
    G_temp=new int[graph_size];
    G_label=new int[graph_size];
    G_order=new int[graph_size];
    G_index=new int[graph_size];
    re_index=new int[graph_size];
    for(int i=0;i<graph_size;++i){
        G_index[i]=-1;
        re_index[i]=-1;
        G_label[i]=0;
        G_temp[i]=degree[i];
        if(degree[i]>max){
            max=degree[i];
        }
    }
    this->max_degree=max;
    
}

void CoreLocate::Coredecompose(){
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
        G_order[node]=i;
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

int CoreLocate::GetMaxcore(){
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
    cout<<"%of max_core: "<<1.0*temp_count/graph_size<<endl;
    cout<<"vertices in graph: "<<graph_size<<endl;
    cout<<"vertices in core: "<<temp_count<<endl;

    return max;
}

int CoreLocate::CorePrune(int **&new_graph, int *&new_degree, int &new_bi){

    int count=0;
    int bi_count=0;
   
    for(int i=0;i<graph_size;++i){
        if(G_label[i]>=K){
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

void CoreLocate::GetReIndex(){
    for(int i=0;i<graph_size;++i){
        if(G_index[i]>=0){
            re_index[G_index[i]]=i;
        }
    }
}

void CoreLocate::CoreOrdering(int **&new_graph, int *&new_degree){
    new_degree=new int[graph_size];
    new_graph=new int*[graph_size];
    for(int i=0;i<graph_size;++i){    
        int *neg=new int[degree[i]];
        new_degree[G_order[i]]=degree[i];
        for(int j=0;j<degree[i];++j){   
            neg[j]=G_order[Graph[i][j]];  
        }
        sort(neg,neg+degree[i]);
        new_graph[G_order[i]]=neg;
    }
}

void CoreLocate::Bipartite_CoreOrdering(int **&new_graph, int *&new_degree){
    new_degree=new int[graph_size];
    new_graph=new int*[graph_size];
    int *pos_index=new int[graph_size];
    for(int i=0;i<bipartite;++i){
        pos_index[i]=0;
    }
    for(int i=bipartite;i<graph_size;++i){
        pos_index[i]=bipartite;
    }

    for(int i=0;i<bipartite;++i){
        for(int j=i+1;j<bipartite;++j){
            if(G_order[i]>G_order[j]){
                pos_index[i]++;
            }else{
                pos_index[j]++;
            }
        }
    }

    for(int i=bipartite;i<graph_size;++i){
        for(int j=i+1;j<graph_size;++j){
            if(G_order[i]>G_order[j]){
                pos_index[i]++;
            }else{
                pos_index[j]++;
            }
        }
    }

    for(int i=0;i<graph_size;++i){    
        int *neg=new int[degree[i]];
        new_degree[pos_index[i]]=degree[i];
        for(int j=0;j<degree[i];++j){   
            neg[j]=pos_index[Graph[i][j]];  
        }
        sort(neg,neg+degree[i]);
        new_graph[pos_index[i]]=neg;
    }
}

void CoreLocate::Bipartite_R_CoreOrdering(int **&new_graph, int *&new_degree){
    new_degree=new int[graph_size];
    new_graph=new int*[graph_size];
    int *pos_index=new int[graph_size];
    for(int i=0;i<bipartite;++i){
        pos_index[i]=0;
    }
    for(int i=bipartite;i<graph_size;++i){
        pos_index[i]=bipartite;
    }

    for(int i=0;i<bipartite;++i){
        for(int j=i+1;j<bipartite;++j){
            if(G_order[i]<G_order[j]){
                pos_index[i]++;
            }else{
                pos_index[j]++;
            }
        }
    }

    for(int i=bipartite;i<graph_size;++i){
        for(int j=i+1;j<graph_size;++j){
            if(G_order[i]<G_order[j]){
                pos_index[i]++;
            }else{
                pos_index[j]++;
            }
        }
    }

    for(int i=0;i<graph_size;++i){    
        int *neg=new int[degree[i]];
        new_degree[pos_index[i]]=degree[i];
        for(int j=0;j<degree[i];++j){   
            neg[j]=pos_index[Graph[i][j]];  
        }
        sort(neg,neg+degree[i]);
        new_graph[pos_index[i]]=neg;
    }
}

