#include<time.h>
#include<fstream>
#include <iostream>
#include <utility>
#include "cmdline.h"
#include "core.h"
#include "biplex1.h"

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
                // cout<<neg[i]<<' ';
            }
            // cout<<endl;
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

pair<int,int> GetMaxPair(map<int,int> pairs){
    int x,y;
    int res = 0;
    for(auto &kv:pairs){
        if(res< kv.first*kv.second){
            x = kv.first;
            y = kv.second;
            res = x*y;
        }
    }
    return make_pair(x,y);
}

int main(int argc, char *argv[])
{
    //parameters
    cmdline::parser a;
    a.add<string>("file", 'f', "filename", true, "");
    a.add<int>("k", 'k', "parameter k", false, 1, cmdline::range(0, 65535));
    a.add<int>("s", 's', "minimum size", false, 3, cmdline::range(1, 65535));
    a.parse_check(argc, argv);

    string filepath = a.get<string>("file");
    int k = a.get<int>("k");
    int s = a.get<int>("s");

    //Read graph
    int bi=0;
    int *degree=NULL;
    int **Graph=NULL;
    int graph_size = ReadGraph(filepath,Graph,degree,bi);
    Core core(Graph, degree, graph_size,bi,k);

    int res = 0;
    int **pG=NULL;
    int *pd=NULL;
    int pb=0;
    int pgs=core.XYCorePrune(pG,pd,pb,s,s);
    cout<<" and core size is "<<pgs<<" and core's bipartite is: "<<pb<<endl;

    time_t begin = clock();
    biplexminer_xy bminer_xy(pG,pgs,pb,pd,s,s,k);
    bminer_xy.Miner();
    time_t end = clock();
    res=bminer_xy.res_size;

    cout<<"result size "<<res<<endl;
    cout<<"Total Time: "<<(double)(end-begin)/CLOCKS_PER_SEC<<"s."<<endl;
    return 0;
}
