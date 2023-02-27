#include <iostream>
#include <bits/stdc++.h>
#include <cstring>
#include <time.h>
#include <chrono>
#include "utils.h"
#include "Corepruning.h"
using namespace std;
bool quiet=false;
int ReadGraph(string dataset_path,int **&Graph, int *&degree, int &bipartite,int &EdgeNUM){
    ifstream read;
    read.open(dataset_path);
    string temp;
    read>>temp;
    cout<<temp<<endl;
    int graph_size=stoi(temp);
    if(!quiet) cout<<"read graph size is "<<graph_size<<endl;
    Graph=new int*[graph_size];
    degree=new int[graph_size];
    read>>temp;
    if(!quiet) cout<<temp<<endl;
    int B_index=stoi(temp);
    bipartite=B_index;
    if(!quiet) cout<<"read bipartite is "<<bipartite<<endl;
    read>>temp;
    int E_index=stoi(temp);
    EdgeNUM=E_index;
    if(!quiet) cout<<"read edgenum is "<<EdgeNUM<<endl;
    int index=0;
    int * neg=new int[graph_size];
    char a;
    int temp_count=0;
    bool first=true;
    while(!read.eof()){
        if(first){
            read>>temp;
            int node=stoi(temp);
            if(!quiet) cout<<"now read node id is "<<node<<" and its neighbor is: ";
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
            if(!quiet) cout<<endl;
            continue;
        }
        read>>temp;
        neg[temp_count]=stoi(temp);
        if(!quiet) cout<<neg[temp_count]<<" ";
        temp_count++;

    }
    delete []neg;
    return (int)graph_size;
}
int main()
{
	string filepath="/home/pardon/DataSet/leadership.txt";
    int k=1;
    int alpha=3;
    int beta=3;
	//Read Graph
	int bi=0;
	int ** Graph=NULL;
	int * Degree=NULL;
	int EdgeNum=0;
	cout<<"now read the data graph "<<endl;
	time_t readgraph1=clock();
	cout<<"filepath is "<<filepath<<endl;
	int Graph_size=ReadGraph(filepath,Graph,Degree,bi,EdgeNum);
	time_t readgraph2=clock();
	cout<<"finish reading graph and time cost: "<<(double)(readgraph2-readgraph1)/CLOCKS_PER_SEC<<"s."<<endl;

    CoreLocate core(Graph,Degree,Graph_size,bi,2*k);
    core.Coredecompose();
    int **pG=NULL;
    int *pd=NULL;
    int pb=0;
    int pgs=core.CorePrune(pG,pd,pb);
	//CSR
	cout<<"now use CSR format to store the graph"<<endl;
	nodeT * Index=new nodeT[GRAPH_SIZE+1];
	nodeT * Neighbor=new nodeT[EDGE_SIZE*2];
	nodeT sum=0;
	Index[0]=0;
	time_t CSR1=clock();
	for(nodeT i=0;i<pgs;i++)
	{
		sum+=pd[i];
		Index[i+1]=sum;
		for(nodeT j=0;j<pd[i];j++)
        {
            Neighbor[Index[i]+j]=pG[i][j];
        }
	}
	time_t CSR2=clock();
	cout<<"CSR Finished! CSR Time: "<<(double)(CSR2-CSR1)/ CLOCKS_PER_SEC<<" s."<<endl;
	if(!quiet)
	{
		cout<<"now print the pruned Graph in CSR format:"<<endl;
		cout<<"index array: ";
		for(int i=0;i<pgs;i++)
		{
			cout<<Index[i]<<" ";
		}
		cout<<endl;
		cout<<"neighbor array: ";
		for(int i=0;i<EdgeNum*2;i++)
		{
			cout<<Neighbor[i]<<" ";
		}
		cout<<endl;
	};
	degreeT * degree=new degreeT[GRAPH_SIZE];
	for(int i=0;i<pgs;i++)
	{
		degree[i]=pd[i];
	}
	nodeT graph_size=pgs;
	nodeT bipartite=pb;
	branch * ddr=new branch[DDR_SIZE];
	branch * mbp=new branch[1];
	nodeT * maximum=new nodeT[1];
	nodeT theta_l=alpha;
	nodeT theta_r=beta;
	cout<<"then deal the data in FPGA FastBB Framwork"<<endl;
	auto framework_begin=std::chrono::steady_clock().now();
	Framework(ddr,Index,Neighbor,degree,graph_size,bipartite,k,alpha,beta,mbp,maximum);
	auto framework_end=std::chrono::steady_clock().now();
	std::chrono::duration<double> framework_time=framework_end-framework_begin;
	cout<<"the whole search time is "<<framework_time.count()<<"sec."<<endl;
	cout<<"maximum biplex size is "<<maximum[0]<<endl;
	cout<<"L("<<mbp[0].P_L.size+mbp[0].C_L.size<<"):";
	for(nodeT i=0;i<mbp[0].P_L.size;i++)
		cout<<mbp[0].P_L.data[i]<<' ';
	for(nodeT i=0;i<mbp[0].C_L.size;i++)
		cout<<mbp[0].C_L.data[i]<<' ';
	cout<<endl;
	cout<<"R("<<mbp[0].P_R.size+mbp[0].C_R.size<<"):";
	for(nodeT i=0;i<mbp[0].P_R.size;i++)
		cout<<mbp[0].P_R.data[i]<<' ';
	for(nodeT i=0;i<mbp[0].C_R.size;i++)
		cout<<mbp[0].C_R.data[i]<<' ';
	cout<<endl;
	return 0;
}
