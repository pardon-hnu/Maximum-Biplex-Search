#pragma once
#include<iostream>
#include <queue>
#include"FastBB.h"
#include <time.h>
using namespace std;

class Framework
{
private:
    int **Graph;
    int *degree;
    int bi_index;
    int graph_size;
    int k;

    int **pG;
    int *pd;
    int pb, pgs;
    int top_K;
    int theta_l,theta_r;

    int *temp_index;
    int *temp_deg;
    int *re_temp_index;

    int temp_i,temp_j,temp_k, temp_node, temp_node2;
    
    void Record_dec();

public:
    Framework(int **Graph, int *degree, int graph_size, int bi_index, int k, int theta_l, int theta_r, int top_K);
    ~Framework();
    void One_hop(int l, int r);
    void PBIE();
    void test(int l,int r);
    int edges,v_r,v_l;
    subgraph MaxBP;
    priority_queue<subgraph,vector<subgraph>,greater<subgraph>> res_top_k;
};

Framework::Framework(int **Graph, int *degree, int graph_size, int bi_index, int k, int theta_l, int theta_r, int top_K)
{
    this->Graph=Graph;
    this->degree=degree;
    this->graph_size=graph_size;
    this->bi_index=bi_index;
    this->k=k;
    this->theta_l=theta_l;
    this->theta_r=theta_r;
    this->top_K=top_K;
    this->pgs=0;

    temp_index=new int[graph_size];
    temp_deg=new int[graph_size];
    re_temp_index=new int[graph_size];

    for(int i=graph_size-1;i>=0;--i){
        temp_index[i]=0;
        temp_deg[i]=0;
        re_temp_index[i]=0;
    }
    
}

Framework::~Framework()
{
}


void Framework::PBIE(){
    int current=0;
    int maxLdeg=0, maxRdeg=0;
    for(temp_i=bi_index-1;temp_i>=0;--temp_i){
        if(degree[temp_i]>maxLdeg){
            maxLdeg=degree[temp_i];
        }
    }
    for(temp_i=graph_size-1;temp_i>=0;--temp_i){
        if(degree[temp_i]>maxRdeg){
            maxRdeg=degree[temp_i];
        }
    }
    int low_l, low_r;
    low_l=maxLdeg;
    int last_v;
    //int bound=2*k+1;
    priority_queue<subgraph,vector<subgraph>,greater<subgraph>> top_k;
    while(low_l>theta_l){
        last_v=low_l;
        low_r=current/low_l<theta_l?theta_l:current/low_l;
        low_l=low_l/2<theta_r?theta_r:low_l/2;
		// last_v=615;
        // low_l=62;
		// low_r=3;
		// current=0;
        cout<<"now low bound of L is "<<low_l<<" and low bound of R is "<<low_r<<endl;
        time_t temptime1=clock();
        One_hop(low_r-k, low_l-k);
        time_t temptime2=clock();
        cout<<(double)(temptime2-temptime1)/CLOCKS_PER_SEC<<endl;
        cout<<pgs<<' '<<pb<<endl;
        if(pgs!=0){
            BKMB bkmb(pG,pd,pgs,pb,k,low_l,low_r,top_K,top_k);
            bkmb.myset=last_v;
            bkmb.edges=current;
            bkmb.re_temp_index=re_temp_index;
            bkmb.BKmb();
            top_k=bkmb.top_k;
            if(!top_k.empty())
            {
                subgraph subg=top_k.top();
                cout<<subg;
                current=bkmb.top_k.top().edges;
            }
            else if(current)
                {current=0;}
        }
        // break;
    }    
    res_top_k=top_k;
}

void Framework::Record_dec(){

    int maxLdeg=0;
    int maxRdeg=0;
    for(temp_i=bi_index-1;temp_i>=0;--temp_i){
        if(degree[temp_i]>maxLdeg){
            maxLdeg=degree[temp_i];
        }
    }
    for(temp_i=graph_size-1;temp_i>=bi_index;--temp_i){
        if(degree[temp_i]>maxRdeg){
            maxRdeg=degree[temp_i];
        }
    }

    int *binL=new int[maxLdeg+1];
    int *binR=new int[maxRdeg+1];

    for(temp_i=maxLdeg;temp_i>=0;--temp_i){
        binL[temp_i]=0;
    }
    for(temp_i=maxRdeg;temp_i>=0;--temp_i){
        binR[temp_i]=0;
    }
    for(temp_i=graph_size-1;temp_i>=bi_index;--temp_i){
        binR[degree[temp_i]]+=1;
    }
    for(temp_i=bi_index-1;temp_i>=0;--temp_i){
        binL[degree[temp_i]]+=1;
    }

    temp_i=0;
    for(temp_j=0;temp_j<=maxLdeg;++temp_j){
        temp_k=binL[temp_j];
        binL[temp_j]=temp_i;
        temp_i+=temp_k;
    }
    temp_i=0;
    for(temp_j=0;temp_j<=maxRdeg;++temp_j){
        temp_k=binR[temp_j];
        binR[temp_j]=temp_i;
        temp_i+=temp_k;
    }

    int *pos=new int[graph_size+1];
    int *vert=new int[graph_size+1];
    int *G_temp=new int[graph_size];

    for(temp_i=graph_size;temp_i>=0;--temp_i){
        pos[temp_i]=0;
        vert[temp_i]=0;

    }
    for(temp_i=graph_size-1;temp_i>=0;--temp_i){
        G_temp[temp_i]=degree[temp_i];
    }
    for(temp_i=0;temp_i<bi_index;++temp_i){
        pos[temp_i]=binL[degree[temp_i]];
        vert[pos[temp_i]]=temp_i;
        binL[degree[temp_i]]+=1;
    }
    for(temp_i=maxLdeg;temp_i>=1;--temp_i){
        binL[temp_i]=binL[temp_i-1];
    }
    binL[0]=1;
    
    for(temp_i=0;temp_i<bi_index;++temp_i){
        temp_node=vert[temp_i];
        for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
            temp_node2=Graph[temp_node][temp_j];
            if(G_temp[temp_node2]>G_temp[temp_node]){
                int du=G_temp[temp_node2]; int pu=pos[temp_node2];
                int pw=binL[du]; int w=vert[pw];
                if(temp_node2!=w){
                    pos[temp_node2]=pw; vert[pu]=w;
                    pos[w]=pu; vert[pw]=temp_node2;
                }
                binL[du]+=1;
                G_temp[temp_node2]-=1;
            }
        }
    }
}

void Framework::One_hop(int l, int r){
    for(temp_i=graph_size-1;temp_i>=0;--temp_i){
        temp_index[temp_i]=1;
        temp_deg[temp_i]=degree[temp_i];
    }
    bool flag=true;
    int countL=0, countR=0;
    while(flag){
        
        flag=false;
        for(temp_i=bi_index-1;temp_i>=0;--temp_i){
            if(temp_index[temp_i]==0||temp_deg[temp_i]>=l) continue;
            temp_index[temp_i]=0;
            flag=true;
            countL++;
            for(temp_j=degree[temp_i]-1;temp_j>=0;--temp_j){
                --temp_deg[Graph[temp_i][temp_j]];
            }
        }
        for(temp_i=bi_index;temp_i<graph_size;++temp_i){
            if(temp_index[temp_i]==0||temp_deg[temp_i]>=r) continue;
            flag=true;
            countR++;
            temp_index[temp_i]=0;
            for(temp_j=degree[temp_i]-1;temp_j>=0;--temp_j){
                --temp_deg[Graph[temp_i][temp_j]];
            }
        }
        
    }
    if(pgs!=0){
        for(temp_i=pgs-1;temp_i>=0;--temp_i){
            delete [] pG[temp_i];
        }
    
        delete []pd;
        delete []pG;
    }
    
    
    pgs=graph_size-countL-countR;
    if(pgs==0) return;
    pd=new int[pgs];
    temp_i=1;

    for(temp_j=0;temp_j<graph_size;++temp_j){
        if(temp_index[temp_j]){
            temp_index[temp_j]=temp_i;
            ++temp_i;
        }
    }

    for(temp_j=0;temp_j<graph_size;++temp_j){
        if(temp_index[temp_j]){
            temp_node=temp_index[temp_j]-1;
            re_temp_index[temp_node]=temp_j;
        }
    }

    pb=bi_index-countL;
    pG=new int*[pgs];
    for(temp_i=0;temp_i<graph_size;++temp_i){
        temp_node=temp_index[temp_i]-1;
        if(temp_node>=0){
            pd[temp_node]=temp_deg[temp_i];
            pG[temp_node]=new int[temp_deg[temp_i]];
            temp_k=0;
            for(temp_j=0;temp_j<degree[temp_i];++temp_j){
                temp_node2=Graph[temp_i][temp_j];
                if(temp_index[temp_node2]){
                    pG[temp_node][temp_k]=temp_index[temp_node2]-1;
                    temp_k++;
                }
            }
        }
    }  
    
}
void Framework::test(int l,int r)
{
    priority_queue<subgraph,vector<subgraph>,greater<subgraph>> top_k;
    One_hop(r-k,l-k);
    if(pgs!=0){
        BKMB bkmb(pG,pd,pgs,pb,k,l,r,top_K,top_k);
        bkmb.myset=pgs;
        bkmb.edges=0;
        bkmb.re_temp_index=re_temp_index;
        bkmb.BKmb();
        res_top_k=bkmb.top_k;
    }
}


