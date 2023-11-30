#include<time.h>
#include<fstream>
#include<cstring>
#include <iostream>
#include <utility>
#include "core.h"
#include "cmdline.h"

using namespace std;

int ReadGraph(string dataset_path,int **&Graph, int *&degree, int &bipartite){
    ifstream read;
    read.open(dataset_path);
    string temp;
    read>>temp;
    int graph_size=stoi(temp);
    Graph=new int*[graph_size];
    delete []degree;
    degree=new int[graph_size];
    read>>temp;
    int B_index=stoi(temp);
    bipartite=B_index;
    read>>temp;
    int index=0;
    int *neg=new int[graph_size];
    char a;
    int temp_count=0;
    bool first=true;
    while(!read.eof()){
        if(first){
            read>>temp;
            first=false;
        }
        read.get(a);
        if(a=='\n'){
            if(index>=graph_size)
                break;
            degree[index]=temp_count;
            int *temp_array=new int[temp_count];
            for(int i=0;i<temp_count;++i){
                temp_array[i]=neg[i];
            }
            Graph[index]=temp_array;
            temp_count=0;
            index++;
            first=true;
            continue;
        }
        read>>temp;
        neg[temp_count]=stoi(temp);
        temp_count++;

    }
    delete []neg;
    return graph_size;
}


int greedy(int **Graph, int *degree, int graph_size, int bipartite, int K){
    int G_temp[graph_size];
    int G_state[graph_size];
    int G_index[graph_size];
    int U = bipartite;
    int V = graph_size-bipartite;

    cout<<"u: "<<U<<", V: "<<V<<", bi: "<<bipartite<<", graph size:"<<graph_size<<endl;

    //重置degree和节点状态
    for(int i=0;i<graph_size;++i){
        G_temp[i]=degree[i];
        G_state[i] = 1;
        G_index[i] = 0;
    }

    for(int i=bipartite;i<graph_size;++i){
        G_index[i] = 1;
    }




    while (1){

        int flag1 = 1;
        int flag2 = 1;


        int du_min = 1000000000;
        int u_min;

        //获取degree最小的节点
        for (int i = 0; i < graph_size; i++){
            if(G_state[i]){
                if(du_min > G_temp[i]){
                    du_min = G_temp[i];
                    u_min = i;
                }
            }
        }

        //remove u_min
        G_state[u_min] = 0;
        for(int i=0;i<degree[u_min];i++){
            G_temp[Graph[u_min][i]] -= 1;
        }
        if(G_index[u_min]){//V
            V--;
        }else{//U
            U--;
        }

        //check
        for(int i=0;i<bipartite;i++){
            if(G_state[i]){
                if(G_temp[i]<V-K)
                    flag1 = 0;
            }
        }
        for(int i=bipartite;i<graph_size;i++){
            if(G_state[i]){
                if(G_temp[i]<U-K)
                    flag2 = 0;
            }
        }


        //cout<<"u: "<<U<<", V: "<<V<<endl;


        if (flag1&&flag2) break;
    }

    //size
    int res_size = 0;
    for (int i = 0; i < bipartite; i++){
        if (G_state[i]){
            res_size += G_temp[i];
        }
    }
    //int res_size = 0;
    cout<<"|U|: "<<U<<", |V|:"<<V<<"."<<endl;
    cout<<"size: "<<res_size<<endl;

    return res_size;
}

int greedy1(int **Graph, int *degree, int graph_size, int bipartite, int K, double eta, int cn){
    cout<<"size: "<<graph_size<<", bi: "<<bipartite<<endl;
    int *G_temp;
    int *G_state;
    G_temp = new int [graph_size];
    G_state = new int [graph_size];
    int U = bipartite;
    int V = graph_size-bipartite;

    cout<<"u: "<<U<<", V: "<<V<<", bi: "<<bipartite<<", graph size:"<<graph_size<<endl;

    //重置degree和节点状态
    for(int i=0;i<graph_size;++i){
        G_temp[i]=degree[i];
        G_state[i] = 1;
    }

    while (1){

        int flag1 = 1;
        int flag2 = 1;


        double a = U;
        double b = V;

        //cout<<"U: "<<a<<", V: "<<b<<endl;

        if(a/b>eta){//移除U中节点
            int du_min = INT32_MAX;
            //int u_min;
            for (int i = 0; i < bipartite; i++){
                if(G_state[i]){
                    if(du_min > G_temp[i]){
                        du_min = G_temp[i];
                        //u_min = i;
                    }
                }
            }

            for (int i = 0; i < bipartite; i++){
                if(G_state[i]){
                    if(du_min >= G_temp[i]){
                        //du_min = G_temp[i];
                        int u_min = i;
                        //remove u_min
                        G_state[u_min] = 0;
                        for(int i=0;i<degree[u_min];i++){
                            G_temp[Graph[u_min][i]] -= 1;
                        }
                        U--;
                    }
                }
            }

            
        }
        else{//移除V中节点
            int dv_min = INT32_MAX;
            //int v_min;
            for (int i = bipartite; i < graph_size; i++){
                if(G_state[i]){
                    if(dv_min > G_temp[i]){
                        dv_min = G_temp[i];
                        //v_min = i;
                    }
                }
            }

            for (int i = bipartite; i < graph_size; i++){
                if(G_state[i]){
                    if(dv_min >= G_temp[i]){
                        //dv_min = G_temp[i];
                        //v_min = i;
                        int v_min = i;

                        //remove v_min
                        G_state[v_min] = 0;
                        for(int i=0;i<degree[v_min];i++){
                            G_temp[Graph[v_min][i]] -= 1;
                        }
                        V--;
                    }
                }
            }

            
        }


        //check
        for(int i=0;i<bipartite;i++){
            if(G_state[i]){
                if(G_temp[i]<V-K)
                    flag1 = 0;
            }
        }
        for(int i=bipartite;i<graph_size;i++){
            if(G_state[i]){
                if(G_temp[i]<U-K)
                    flag2 = 0;
            }
        }


        //cout<<"u: "<<U<<", V: "<<V<<endl;


        if (flag1&&flag2) break;
    }

    //size
    int res_size = 0;
    for (int i = 0; i < bipartite; i++){
        if (G_state[i]){
            res_size += G_temp[i];
        }
    }
    if(U<cn||V<cn){
        cout<<"|U|: "<<U<<", |V|:"<<V<<"."<<endl;
        cout<<"size: "<<res_size;
        cout<<", no result."<<endl;
        return 0;
    }else{
        cout<<"|U|: "<<U<<", |V|:"<<V<<"."<<endl;
        cout<<"U: ";
        for(int i=0;i<bipartite;i++)
        {
            if(G_state[i])
            cout<<i<<":"<<degree[i]<<" "" ";
        }
        cout<<endl;
        cout<<"V: ";
        for(int i=bipartite;i<graph_size;i++)
        {
            if(G_state[i])
            cout<<i<<":"<<degree[i]<<" ";
        }
        cout<<endl; 
        cout<<"size: "<<res_size<<endl;
    }
    return res_size;
}


int main(int argc, char *argv[])
{
    //parameters
    cmdline::parser a;
    a.add<string>("file", 'f', "filename", true, "");
    a.add<int>("k", 'k', "parameter k", false, 1, cmdline::range(1, 65535));
    a.add<int>("s", 's', "minimum size", false, 7, cmdline::range(1, 65535));
    a.parse_check(argc, argv);

    string filepath = a.get<string>("file");
    int k = a.get<int>("k");
    int s = a.get<int>("s");

    //Read graph
    int bi=0;
    int *degree=NULL;
    int **Graph=NULL;
    int graph_size = ReadGraph(filepath,Graph,degree,bi);

    time_t begin = clock();
    Core core(Graph, degree, graph_size,bi,k);
    core.Coredecompose();
    int maxcore = core.GetMaxcore();
    int res = 0;

    int cn = 7;
    cout<<"max core: "<<maxcore<<endl;
    int **pG=NULL;
    int *pd=NULL;
    int pb=0;
    int pgs=core.CorePrune(pG,pd,pb,cn-k);
    cn=s;
    double eta = 10;
    while(eta){
        int temp_res = greedy1(pG,pd,pgs,pb,k,eta,cn);
        if(res < temp_res)
            res = temp_res;
        if(eta>1)
            eta -= 0.5;
        else
            eta -= 0.05;
        if(eta < 0)
            break;
    }
    time_t end = clock();
    
    cout<<"Total Time: "<<(double)(end-begin)/ CLOCKS_PER_SEC<<" s."<<endl;
    cout<<"Max size: "<<res<<endl;
    return 0;
}
