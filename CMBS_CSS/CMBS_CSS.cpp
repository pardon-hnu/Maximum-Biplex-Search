#include<time.h>
#include<fstream>
#include <iostream>
#include <utility>
#include "cmdline.h"
#include "core.h"
#include "FastBB.h"

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
    Core core(Graph, degree, graph_size,bi,k);

    //int x = 4;
    //int y = core.GetMaxY(x);
    //printf("x: %d, maxy: %d.\n",x,y);

        core.Coredecompose();
        int maxcore = core.GetMaxcore();
        int cn = maxcore;
        int res = 0;
        int res_size = 0;
        //printf("max core number: %d.\n", maxcore);
        //printf("size: %d.\n",pgs);

        //biplex calculation
        time_t begin = clock();
        int **pG=NULL;
        int *pd=NULL;
        int pb=0;
        //compute xy-cores from s to cn
        for(int x = core.GetMaxX(s);x>=s;x--){
            //initial y*
            cout<<"now x is "<<x<<endl;
            int y = core.GetMaxY(x);
            int bound = res_size/x;
            //compute xy-core contain biplex
            while(1){
                if(y<=bound||y<s) break;
                //subgraph of xy-core
                int **pG=NULL;
                int *pd=NULL;
                int pb=0;
                //xy-core purning
                time_t s1 = clock();
                int pgs=core.XYCorePrune(pG,pd,pb,y-k,x-k);
                time_t s2 = clock();

                cout<<"x: "<<x<<" y: "<<y<<" and core size is "<<pgs<<" and core's bipartite is: "<<pb<<endl;

                //compute biplex
                time_t s3 = clock();
                BKMB bkmb(pG,pd,pgs,pb,k,x,y);
                bkmb.edges=s*(s-k);
                bkmb.BKmb2();
        
                //query
                if(bkmb.edges>res_size){
                    res_size = bkmb.edges;
                    cout<<"size: "<<res_size<<endl;
                    break;
                }
                else
                    y -= 1;
            
                time_t s4 = clock();
                cout<<"core Time: "<<(double)(s2-s1)/CLOCKS_PER_SEC<<". biplex Time: "<<(double)(s4-s3)/CLOCKS_PER_SEC<<"s."<<endl;

            }
        }
    time_t end2 = clock();
    cout<<"Total Time: "<<(double)(end2-begin)/ CLOCKS_PER_SEC<<" s."<<endl;
    cout<<"Size: "<<res_size<<endl;
    return 0;
}
