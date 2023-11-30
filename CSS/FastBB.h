#pragma once
#include<iostream>
#include<vector>
#include<time.h>
#include"RandList.h"
#include"subgraph.h"
#include<queue>
using namespace std;
#define DEBUG 0
class BKMB
{
private:
    /* data */
    int **Graph;
    int *degree;
    int graph_size;
    int bi_index;
    int k;
    int left, right;
    int temp_i,temp_j,temp_k, temp_value, temp_node, temp_node2, temp_size1, temp_size2;
    MBitSet *bit_G;
    RandList P_L, P_R, Cand_L, Cand_R, Exc_L, Exc_R;
    int *degInP;
    int *degInG;
    int *G_index; int count;
    vector<int> temp_vector;
    int topK;
    int bs;
    bool CanBranch();
    void BKmb_Rec();
    void Record(vector<int>& left, vector<int>& right);
    void Record_now(vector<int>& left, vector<int>& right);
public:
    BKMB(int **Graph, int *degree, int graph_size, int bi_index, int k,int left, int right, int topK, priority_queue<subgraph,vector<subgraph>,greater<subgraph>> top_k);
    ~BKMB();
    void BKmb();
    int total_num;
    int tt;
    int edges,v_l,v_r;
    int myset;
    //#ifdef OUTPUT
    int *re_temp_index;
    // subgraph *maxbp;
    priority_queue<subgraph,vector<subgraph>,greater<subgraph>> top_k;
    //#endif
    bool debug;
};

BKMB::BKMB(int **Graph, int *degree, int graph_size, int bi_index, int k, int left, int right, int topK, priority_queue<subgraph,vector<subgraph>,greater<subgraph>> top_k)
{
    this->Graph=Graph;
    this->degree=degree;
    this->graph_size=graph_size;
    this->bi_index=bi_index;
    this->k=k;
    this->left=left;
    this->right=right;
    total_num=0;
    edges=0;
    myset=graph_size;
    bs=0;
    P_L.init(graph_size);
    P_R.init(graph_size);
    Cand_L.init(graph_size);
    Cand_R.init(graph_size);
    Exc_L.init(graph_size);
    Exc_R.init(graph_size);
    temp_vector.reserve(graph_size);
    tt=0;
    this->top_k=top_k;
    this->topK=topK;
    if(top_k.empty()) edges=0;
    else{
        subgraph temp=top_k.top();
        edges=temp.edges;
    }

    degInG=new int[graph_size];
    degInP=new int[graph_size];
    G_index=new int[graph_size];
    count=0;

    for(temp_i=graph_size-1;temp_i>=0;--temp_i){
        degInP[temp_i]=0;
        degInG[temp_i]=0;
        G_index[temp_i]=0;
    }
    

    bit_G=new MBitSet[graph_size];
    double ttm=0;
    for(int i=0;i<graph_size;++i){
        bit_G[i].allocacte(graph_size+1);
        bit_G[i].reinit(graph_size+1);
        for(int j=0;j<degree[i];++j){
            bit_G[i].set(Graph[i][j]);
        }
        ttm+=sizeof(bit_G[0]);
    }
}

BKMB::~BKMB()
{
    delete [] degInG;
    delete [] degInP;
    delete [] G_index;
    //delete []bit_G;
    //delete &temp_vector;
    //delete &P_L, &P_R, &Cand_L, &Cand_R, &Exc_L, &Exc_R;
}

void BKMB::Record(vector<int>& left, vector<int>& right){
    for(int i=P_L.vnum-1;i>=0;--i)
        left.push_back(re_temp_index[P_L.vlist[i]]);
    for(int i=Cand_L.vnum-1;i>=0;--i)
        left.push_back(re_temp_index[Cand_L.vlist[i]]);
    for(int i=P_R.vnum-1;i>=0;--i)
        right.push_back(re_temp_index[P_R.vlist[i]]);
    for(int i=Cand_R.vnum-1;i>=0;--i)
        right.push_back(re_temp_index[Cand_R.vlist[i]]);
}
void BKMB::Record_now(vector<int>& left, vector<int>& right){
    for(int i=P_L.vnum-1;i>=0;--i)
        left.push_back(P_L.vlist[i]);
    for(int i=Cand_L.vnum-1;i>=0;--i)
        left.push_back(Cand_L.vlist[i]);
    for(int i=P_R.vnum-1;i>=0;--i)
        right.push_back(P_R.vlist[i]);
    for(int i=Cand_R.vnum-1;i>=0;--i)
        right.push_back(Cand_R.vlist[i]);
}
void BKMB::BKmb(){
    bool flag;
    int sum=0;
    for(int i=0;i<bi_index;++i){
    // each node i in L to select B_i
        // i=3;
        P_L.add(i);
        for(temp_i=degree[i]-1;temp_i>=0;--temp_i){
            temp_node=Graph[i][temp_i]; //u in N(v)
            count=0;
            for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                temp_node2=Graph[temp_node][temp_j];    //v' in N(N(v))
                if(temp_node2>=i) count++;              //v' future in Cand_L must > v (exclude v'<v)
            }
            if(count<left-k) continue;                  //if d(u,Cand_L and v)<left-k  u is prunned
            degInG[temp_node]++;                        
            degInP[temp_node]++;
            Cand_R.add(temp_node);                      //if satisfy, add u to Cand_R and d(u) has connect to v so d(u)++
            for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                temp_node2=Graph[temp_node][temp_j];    
                degInG[temp_node2]++;                   //modify the degree of v'
            }
        }

        //branch and bound && theta-2k constarint
        for(temp_i=degree[i]-1;temp_i>=0;--temp_i){
            temp_node=Graph[i][temp_i];                 //u in N(v)
            for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                temp_node2=Graph[temp_node][temp_j];    //v' in N(N(v))
                if(temp_node2>i&&degInG[temp_node2]>=right-2*k&&Cand_L.vpos[temp_node2]==graph_size){
                /* v' must satisfy three condition:
                    1. v' > v (branch and bound)
                    2. d(v',N(v))>= theta-2k
                    3. v' is not in Cand_L
                */
                    Cand_L.add(temp_node2);             //add Cand_L and update the degree
                    for(temp_k=degree[temp_node2]-1;temp_k>=0;--temp_k){
                        degInG[Graph[temp_node2][temp_k]]++;
                    }
                }
                if(temp_node2<i&&degInG[temp_node2]>=right-2*k&&Exc_L.vpos[temp_node2]==graph_size){
                    Exc_L.add(temp_node2);
                }
                // v' is satisfy the theta-2k constaint but < v added into Exc_l
            }
        }
        // L-side is ok then deal R-side
        for(temp_i=bi_index;temp_i<graph_size;++temp_i){
            if(degInG[temp_i]>=left-k&&Cand_R.vpos[temp_i]==graph_size){
                // if u not in Cand_R, d(u,Cand_L and v)>=left-k
                // add u into Cand_R and update the degree 
                Cand_R.add(temp_i);
                for(temp_j=degree[temp_i]-1;temp_j>=0;--temp_j){
                    degInG[Graph[temp_i][temp_j]]++;
                }                
            }
            if(degInG[temp_i]<left-k&&Cand_R.vpos[temp_i]<graph_size){
                // if u in Cand_R but d(u,Cand_L and v)<left-k
                //remove the u and update the degree
                Cand_R.remove(temp_i);
                for(temp_j=degree[temp_i]-1;temp_j>=0;--temp_j){
                    temp_node=Graph[temp_i][temp_j];
                    degInG[temp_node]--;
                    //update the degree of v' in Cand_L 
                    //if degree<right-k remove and update the degree
                    if(degInG[temp_node]<right-k&&Cand_L.vpos[temp_node]<graph_size){
                        Cand_L.remove(temp_node);
                        for(temp_k=degree[temp_node]-1;temp_k>=0;--temp_k){
                            --degInG[Graph[temp_node][temp_k]];
                        }
                    }
                }
            }
        }
        //pruning 2
        for(temp_i=bi_index;temp_i<graph_size;++temp_i){
            if(degInG[temp_i]>=left-k&&Cand_R.vpos[temp_i]==graph_size){
                Cand_R.add(temp_i);
                for(temp_j=degree[temp_i]-1;temp_j>=0;--temp_j){
                    degInG[Graph[temp_i][temp_j]]++;
                }
            }
            if(degInG[temp_i]<left-k&&Cand_R.vpos[temp_i]<graph_size){
                Cand_R.remove(temp_i);
                for(temp_j=degree[temp_i]-1;temp_j>=0;--temp_j){
                    temp_node=Graph[temp_i][temp_j];
                    degInG[temp_node]--;
                    if(degInG[temp_node]<right-k&&Cand_L.vpos[temp_node]<graph_size){
                        Cand_L.remove(temp_node);
                        for(temp_k=degree[temp_node]-1;temp_k>=0;--temp_k){
                            --degInG[Graph[temp_node][temp_k]];
                        }
                    }
                }
            }
        }
        //pruning 3
        for(temp_i=bi_index;temp_i<graph_size;++temp_i){
            if(degInG[temp_i]>=left-k&&Cand_R.vpos[temp_i]==graph_size){
                Cand_R.add(temp_i);
                for(temp_j=degree[temp_i]-1;temp_j>=0;--temp_j){
                    degInG[Graph[temp_i][temp_j]]++;
                }
            }
            if(degInG[temp_i]<left-k&&Cand_R.vpos[temp_i]<graph_size){
                Cand_R.remove(temp_i);
                for(temp_j=degree[temp_i]-1;temp_j>=0;--temp_j){
                    temp_node=Graph[temp_i][temp_j];
                    degInG[temp_node]--;
                    if(degInG[temp_node]<right-k&&Cand_L.vpos[temp_node]<graph_size){
                        Cand_L.remove(temp_node);
                        for(temp_k=degree[temp_node]-1;temp_k>=0;--temp_k){
                            --degInG[Graph[temp_node][temp_k]];
                        }
                    }
                }
            }
        }


        // v' in exclude set: if v' connect all nodes in Cand_R,the branch is not maximum
        flag=true;
        for(temp_i=Exc_L.vnum-1;temp_i>=0;--temp_i){
            temp_node=Exc_L.vlist[temp_i];
            if(degInG[temp_node]==Cand_R.vnum){
                flag=false;
                break;
            }
        }
        

        if(flag)
        {if(DEBUG) {cout<<"now maximum is"<<edges<<endl;cout<<"("<<i<<")"<<"b0 can be expand"<<endl;}BKmb_Rec();}

        for(temp_i=Cand_R.vnum-1;temp_i>=0;--temp_i){
            temp_node=Cand_R.vlist[temp_i];
            for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                --degInG[Graph[temp_node][temp_j]];
            }
        }
        for(temp_i=Cand_L.vnum-1;temp_i>=0;--temp_i){
            temp_node=Cand_L.vlist[temp_i];
            for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                --degInG[Graph[temp_node][temp_j]];
            }
        }
        
        for(temp_i=P_L.vnum-1;temp_i>=0;--temp_i){
            temp_node=P_L.vlist[temp_i];
            for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                --degInG[Graph[temp_node][temp_j]];
                --degInP[Graph[temp_node][temp_j]];
            }
        }
        
        Cand_R.clear();
        Cand_L.clear();
        P_L.clear();
        Exc_L.clear();   
        // break;
    }
    cout<<"branch num is "<<bs<<endl;
}

void BKMB::BKmb_Rec(){
    debug=false;
    bs++;
    // if(P_L.vnum==1&&P_R.vnum==1&&Cand_L.vnum==152&&Cand_R.vnum==12) debug=true;
    //myset is the upper bound of the L side 
    if(P_L.vnum>=myset){
        return;
    }// if |L|> upper_l return
    
    if(P_L.vnum+Cand_L.vnum<left||P_R.vnum+Cand_R.vnum<right) return;
    //if |L|<low_l or |R| < low_r
    
    if((P_L.vnum+Cand_L.vnum)*(P_R.vnum+Cand_R.vnum)<=edges&&(P_L.vnum+Cand_L.vnum)*(P_R.vnum+Cand_R.vnum)>0) return;
    // if there are no possibility > current maximum return

    // if there are no possibility > current maximum return
    temp_j=0;
    for(temp_i=P_L.vnum-1;temp_i>=0;--temp_i){
        temp_node=P_L.vlist[temp_i];
        temp_j+=degInG[temp_node];
    }
    for(temp_i=Cand_L.vnum-1;temp_i>=0;--temp_i){
        temp_node=Cand_L.vlist[temp_i];
        temp_j+=degInG[temp_node];
    }
    if(temp_j<=edges){
        return;
    }
    
    // if there don't have nodes in Candidate sets
    if(Cand_R.empty()&&Cand_L.empty()){
        if(Exc_R.empty()&&Exc_L.empty()){
            total_num++;
            temp_j=0;
            for(temp_i=P_L.vnum-1;temp_i>=0;--temp_i){
                temp_node=P_L.vlist[temp_i];
                temp_j+=degInP[temp_node];
            }
            //?when now edges > current edges
            if(edges<temp_j){
                if(DEBUG){
                cout<<"find a biplex and size is "<<temp_value<<" and now maximum is "<<edges<<endl;
                cout<<"P_L("<<P_L.vnum<<"):";
                for(int i=0;i<P_L.vnum;i++)
                    cout<<P_L.vlist[i]<<' ';
                cout<<endl;
                cout<<"P_R("<<P_R.vnum<<"):";
                for(int i=0;i<P_R.vnum;i++)
                    cout<<P_R.vlist[i]<<' ';
                cout<<endl;
                cout<<"C_L("<<Cand_L.vnum<<"):";
                for(int i=0;i<Cand_L.vnum;i++)
                    cout<<Cand_L.vlist[i]<<' ';
                cout<<endl;
                cout<<"C_R("<<Cand_R.vnum<<"):";
                for(int i=0;i<Cand_R.vnum;i++)
                    cout<<Cand_R.vlist[i]<<' ';
                cout<<endl;}
                if(top_k.size()<topK){ //if(there has <K results) add now result to top_k
                    subgraph temp_sub(temp_j,P_L.vnum+Cand_L.vnum,P_R.vnum+Cand_R.vnum);
                    Record(temp_sub.left,temp_sub.right);
                    top_k.push(temp_sub);
                }else if(edges<temp_j){ //else if add now result to top_k top_k delete the miner one update the edges
                    subgraph temp_sub(temp_j,P_L.vnum+Cand_L.vnum,P_R.vnum+Cand_R.vnum);
                    Record(temp_sub.left,temp_sub.right);
                    top_k.push(temp_sub);
                    top_k.pop();
                    edges=top_k.top().edges;
                }
            }
        } 
        return;
    }
    if(DEBUG)
    {
        RandList te;te.init(graph_size);
        for(int i=0;i<P_L.vnum;i++)
        {
            te.add(P_L.vlist[i]);
        }
        sort(te.vlist,te.vlist+te.vnum);
        cout<<"P_L("<<te.vnum<<"):";
        for(int i=0;i<te.vnum;i++)
            cout<<te.vlist[i]<<' ';
        cout<<endl;
        te.init(graph_size);
        for(int i=0;i<P_R.vnum;i++)
        {
            te.add(P_R.vlist[i]);
        }
        sort(te.vlist,te.vlist+te.vnum);
        cout<<"P_R("<<te.vnum<<"):";
        for(int i=0;i<te.vnum;i++)
            cout<<te.vlist[i]<<' ';
        cout<<endl;
        te.init(graph_size);
        for(int i=0;i<Cand_L.vnum;i++)
        {
            te.add(Cand_L.vlist[i]);
        }
        sort(te.vlist,te.vlist+te.vnum);
        cout<<"C_L("<<Cand_L.vnum<<"):";
        for(int i=0;i<te.vnum;i++)
            cout<<te.vlist[i]<<' ';
        cout<<endl;
        te.init(graph_size);
        for(int i=0;i<Cand_R.vnum;i++)
        {
            te.add(Cand_R.vlist[i]);
        }
        sort(te.vlist,te.vlist+te.vnum);
        cout<<"C_R("<<Cand_R.vnum<<"):";
        for(int i=0;i<te.vnum;i++)
            cout<<te.vlist[i]<<' ';
        cout<<endl;
    }
    int povit=-1;
    temp_value=0;
    temp_node2=graph_size;
    temp_j=P_L.vnum+Cand_L.vnum;
    temp_k=P_R.vnum+Cand_R.vnum;

    for(temp_i=P_L.vnum-1;temp_i>=0;--temp_i){
        temp_node=P_L.vlist[temp_i];
        if(temp_k-degInG[temp_node]>temp_value){
            temp_value=temp_k-degInG[temp_node];
            povit=temp_node;
        }
        if(temp_node2>degInG[temp_node]){
            temp_node2=degInG[temp_node];
        }
    }// temp_node2:the d_min(v) in L, temp_value:the maximum disconnect to v in L  povit: node v
    if(temp_node2+k<right){
        return;
    }//if the d_min(v)+k<low_r, the branch is prunned
    temp_size1=graph_size;
    for(temp_i=P_R.vnum-1;temp_i>=0;--temp_i){
        temp_node=P_R.vlist[temp_i];
        if(temp_j-degInG[temp_node]>temp_value){
            temp_value=temp_j-degInG[temp_node];
            povit=temp_node;
        }
        if(temp_size1>degInG[temp_node]){
            temp_size1=degInG[temp_node];
        }
    }//temp_size1:the d_min in R, update the temp_value and provit
    if(temp_size1+k<left||(temp_node2+k)*(temp_size1+k)<=edges+2*k){
        return;
    }//if the d_min(u)+k<low_l or upper bound of this banch lower than the current edges,the branch is prunned
    
    if(povit>=0&&temp_value>k){
        if(DEBUG) cout<<"find a pivot in P : "<<povit<<endl;
        RandList *Cand1, *Cand2, *P1, *P2, *Exc1, *Exc2;
        int temp_stand1, temp_stand2;
        if(povit<bi_index){
            Cand1=&Cand_L;
            Cand2=&Cand_R;
            P1=&P_L;
            P2=&P_R;
            Exc1=&Exc_L;
            Exc2=&Exc_R;
            temp_stand1=left;
            temp_stand2=right;
        }else{
            Cand1=&Cand_R;
            Cand2=&Cand_L;
            P1=&P_R;
            P2=&P_L;
            Exc1=&Exc_R;
            Exc2=&Exc_L;
            temp_stand1=right;
            temp_stand2=left;
        }

        vector<int> doing;
        doing.reserve(Cand2->vnum);
        for(temp_i=Cand2->vnum-1;temp_i>=0;--temp_i){
            temp_node=Cand2->vlist[temp_i];
            if(!bit_G[povit].test(temp_node)){
                doing.push_back(temp_node);
            }
        }// find the node disconnect to the povit : gamma(povit,cand)
        int p=k-(P2->vnum-degInP[povit]), idx=1;
        //  the max number nodes in gamma(povit,cand) can be added
        int rec_a, rec_b;
        rec_a=doing[0];
        Cand2->remove(rec_a);
        Exc2->add(rec_a);
        for(temp_i=degree[rec_a]-1;temp_i>=0;--temp_i){
            --degInG[Graph[rec_a][temp_i]];
        }
        //first branch  s+=NULL,c remove first element,d add first element
        if(CanBranch()&&P2->vnum+Cand2->vnum>=temp_stand2)
            {if(!debug)BKmb_Rec();} //if canbranch: there is no node n in D,make branch not maximum
        
        Exc2->remove(rec_a);
        vector<int> remove_C1, remove_E1, remove_C2, remove_E2;
        remove_C1.reserve(Cand1->vnum);
        remove_C2.reserve(Cand2->vnum);
        remove_E1.reserve(Exc1->vnum);
        remove_E2.reserve(Exc2->vnum);
      
        //while add node which in gamma(povit,cand) to branch don't violate the k-biplex definition
        while(idx<p+1){
            P2->add(rec_a);
            for(temp_i=degree[rec_a]-1;temp_i>=0;--temp_i){
                temp_node=Graph[rec_a][temp_i];
                degInG[temp_node]++;
                degInP[temp_node]++;
            }

            if(degInP[rec_a]==P1->vnum-k){
                for(temp_i=Cand1->vnum-1;temp_i>=0;--temp_i){
                    temp_node=Cand1->vlist[temp_i];
                    if(!bit_G[rec_a].test(temp_node)||degInG[temp_node]<temp_stand2-k){
                        Cand1->remove(temp_node);
                        remove_C1.push_back(temp_node);
                        for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                            temp_node2=Graph[temp_node][temp_j];
                            --degInG[temp_node2];
                        }
                    }
                    
                }

                for(temp_i=Exc1->vnum-1;temp_i>=0;--temp_i){
                    temp_node=Exc1->vlist[temp_i];
                    if(!bit_G[rec_a].test(temp_node)||degInG[temp_node]<temp_stand2-k){
                        Exc1->remove(temp_node);
                        remove_E1.push_back(temp_node);
                    }
                }
            }else{
                temp_size1=P2->vnum-k;
                for(temp_i=Cand1->vnum-1;temp_i>=0;--temp_i){
                    temp_node=Cand1->vlist[temp_i];
                    if(degInP[temp_node]<temp_size1||degInG[temp_node]<temp_stand2-k){
                        Cand1->remove(temp_node);
                        remove_C1.push_back(temp_node);
                        for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                            temp_node2=Graph[temp_node][temp_j];
                            --degInG[temp_node2];
                        }
                    }
                }

                for(temp_i=Exc1->vnum-1;temp_i>=0;--temp_i){
                    temp_node=Exc1->vlist[temp_i];
                    if(degInP[temp_node]<temp_size1||degInG[temp_node]<temp_stand2-k){
                        Exc1->remove(temp_node);
                        remove_E1.push_back(temp_node);
                    }
                }
            }
            
            if(debug)
            {
                cout<<"first judge and pruning is finish"<<endl;
                cout<<"C_L("<<Cand_L.vnum<<"):";
                for(int i=0;i<Cand_L.vnum;i++)
                    cout<<Cand_L.vlist[i]<<' ';
                cout<<endl;
                cout<<"C_R("<<Cand_R.vnum<<"):";
                for(int i=0;i<Cand_R.vnum;i++)
                    cout<<Cand_R.vlist[i]<<' ';
                cout<<endl;
            }

            if(Cand1->vnum+P1->vnum<temp_stand1){
                break;
            }

            temp_vector.clear();
            temp_size1=P2->vnum;
            temp_size2=P1->vnum;
            if(debug) cout<<"temp_size1: "<<temp_size1<<" temp_size2: "<<temp_size2<<endl;
            for(temp_i=P1->vnum-1;temp_i>=0;--temp_i){
                temp_node=P1->vlist[temp_i];
                if(degInP[temp_node]==temp_size1-k&&!bit_G[rec_a].test(temp_node)){
                    temp_vector.push_back(temp_node);
                }
            }
            for(temp_i=degree[rec_a]-1;temp_i>=0;--temp_i){
                temp_node=Graph[rec_a][temp_i];
                if(Cand1->contains(temp_node)||P1->contains(temp_node))
                    G_index[temp_node]=1;
            }

            if(debug)
            {
                cout<<"temp_vector("<<temp_vector.size()<<"):";
                for(int bi=0;bi<temp_vector.size();bi++)
                    cout<<temp_vector[bi]<<' ';
                cout<<endl;
                cout<<"G_index: ";
                for(int bi=0;bi<graph_size;bi++)
                {
                    if(G_index[bi]==1) cout<<bi<<' ';
                }
                cout<<endl;
                cout<<"degree: "<<endl;
                for(int bi=0;bi<graph_size;bi++)
                {
                    if(degInG[bi]||degInP[bi])
                        cout<<bi<<":("<<degInP[bi]<<","<<degInG[bi]<<")"<<endl;
                }
            }        
            
            if(!temp_vector.empty()){
                for(temp_i=Cand2->vnum-1;temp_i>=0;--temp_i){
                    temp_node=Cand2->vlist[temp_i];
                    if(debug) cout<<temp_node<<" "<<degInP[temp_node]<<" "<<degInG[temp_node]<<endl;
                    if(degInP[temp_node]==temp_size2) continue;
                    if(degInP[temp_node]<temp_size2-k||degInG[temp_node]<temp_stand1-k){
                        if(debug) cout<<"remove"<<endl;
                        Cand2->remove(temp_node);
                        remove_C2.push_back(temp_node);
                        for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                            temp_node2=Graph[temp_node][temp_j];
                            --degInG[temp_node2];
                        }
                    }else{
                        for(int i:temp_vector){
                            if(!bit_G[temp_node].test(i)){
                                Cand2->remove(temp_node);
                                remove_C2.push_back(temp_node);
                                for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                                    temp_node2=Graph[temp_node][temp_j];
                                    --degInG[temp_node2];
                                }
                                break;
                            }
                        }

                        if(Cand2->contains(temp_node)){
                            count=0;
                            for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                                if(G_index[Graph[temp_node][temp_j]]){
                                    count++;
                                }
                            }
                            if(count<temp_stand1-2*k){
                                Cand2->remove(temp_node);
                                remove_C2.push_back(temp_node);
                                for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                                    temp_node2=Graph[temp_node][temp_j];
                                    --degInG[temp_node2];
                                }
                            }
                        }   
                    }
                }

                for(temp_i=Exc2->vnum-1;temp_i>=0;--temp_i){
                    temp_node=Exc2->vlist[temp_i];
                    if(degInP[temp_node]==temp_size2) continue;
                    if(degInP[temp_node]<temp_size2-k||degInG[temp_node]<temp_stand1-k){
                        Exc2->remove(temp_node);
                        remove_E2.push_back(temp_node);
                    }else{
                        for(int i:temp_vector){
                            if(!bit_G[temp_node].test(i)){
                                Exc2->remove(temp_node);
                                remove_E2.push_back(temp_node);
                                break;
                            }
                        }
                    }
                }
            }

            if(debug)
            {
                cout<<"second judge and pruning is finish"<<endl;
                cout<<"C_L("<<Cand_L.vnum<<"):";
                for(int i=0;i<Cand_L.vnum;i++)
                    cout<<Cand_L.vlist[i]<<' ';
                cout<<endl;
                cout<<"C_R("<<Cand_R.vnum<<"):";
                for(int i=0;i<Cand_R.vnum;i++)
                    cout<<Cand_R.vlist[i]<<' ';
                cout<<endl;
            }

            for(temp_i=degree[rec_a]-1;temp_i>=0;--temp_i){
                temp_node=Graph[rec_a][temp_i];
                G_index[temp_node]=0;
            }
            

            if(Cand2->vnum+P2->vnum<temp_stand2){
                break;
            }

            for(;idx<doing.size()&&!Cand2->contains(doing[idx]);++p,++idx);
            if(idx==doing.size()) break;
            rec_b=doing[idx];

            Cand2->remove(rec_b);
            Exc2->add(rec_b);
            for(temp_i=degree[rec_b]-1;temp_i>=0;--temp_i){
                temp_node=Graph[rec_b][temp_i];
                --degInG[temp_node];
            }
            if(CanBranch())
            {
            if(debug){
                cout<<"*********"<<endl;
            cout<<"P_L("<<P_L.vnum<<"):";
            for(int i=0;i<P_L.vnum;i++)
                cout<<P_L.vlist[i]<<' ';
            cout<<endl;
            cout<<"P_R("<<P_R.vnum<<"):";
            for(int i=0;i<P_R.vnum;i++)
                cout<<P_R.vlist[i]<<' ';
            cout<<endl;
            cout<<"C_L("<<Cand_L.vnum<<"):";
            for(int i=0;i<Cand_L.vnum;i++)
                cout<<Cand_L.vlist[i]<<' ';
            cout<<endl;
            cout<<"C_R("<<Cand_R.vnum<<"):";
            for(int i=0;i<Cand_R.vnum;i++)
                cout<<Cand_R.vlist[i]<<' ';
            cout<<endl;
            cout<<"*********"<<endl;}
            if(!debug)BKmb_Rec();}

            Exc2->remove(rec_b);
            rec_a=rec_b;
            ++idx;
        }
        if(debug) cout<<Cand2->vnum<<endl;
        if(debug) cout<<"idx1: "<<idx<<endl;
        if(Cand1->vnum+P1->vnum>=temp_stand1&&Cand2->vnum+P2->vnum>=temp_stand2){
            
            count=0;
            for(;idx<doing.size();++idx){
                temp_node=doing[idx];
                if(Cand2->contains(temp_node)){
                    count++;
                }
            }
            if(debug) cout<<"idx2: "<<idx<<endl;
            if(Cand2->vnum+P2->vnum-count>=temp_stand2){
                for(;idx<doing.size();++idx){
                temp_node=doing[idx];
                if(Cand2->contains(temp_node)){
                    Cand2->remove(temp_node);
                    for(temp_i=degree[temp_node]-1;temp_i>=0;--temp_i){
                        --degInG[Graph[temp_node][temp_i]];
                    }
                }
                }
                if(debug) cout<<"idx3: "<<idx<<endl;
                if(CanBranch())
                    {BKmb_Rec();}
            }
            
            
            
        }
        


        
        for(int i:doing){
            if(P2->vpos[i]<graph_size){
                P2->remove(i);
                Cand2->add(i);
                for(temp_i=degree[i]-1;temp_i>=0;--temp_i){
                    temp_node=Graph[i][temp_i];
                    --degInP[temp_node];
                }
            }
            if(Cand2->vpos[i]>=graph_size){
                Cand2->add(i);
                for(temp_i=degree[i]-1;temp_i>=0;--temp_i){
                    temp_node=Graph[i][temp_i];
                    degInG[temp_node]++;
                }

            }
        }

        for(int i:remove_C2){
            if(Cand2->vpos[i]==graph_size){
                Cand2->add(i);
                for(temp_i=degree[i]-1;temp_i>=0;--temp_i){
                    temp_node=Graph[i][temp_i];
                    degInG[temp_node]++;
                }
            }
        }

        for(int i:remove_C1){
            Cand1->add(i);
            for(temp_i=degree[i]-1;temp_i>=0;--temp_i){
                temp_node=Graph[i][temp_i];
                degInG[temp_node]++;
            }
        }

        for(int i:remove_E1){
            Exc1->add(i);
        }
        for(int i:remove_E2){
            Exc2->add(i);
        }

        

    }else{
        for(temp_i=Cand_L.vnum-1;temp_i>=0;--temp_i){
            temp_node=Cand_L.vlist[temp_i];
            if(temp_k-degInG[temp_node]>temp_value){
                temp_value=temp_k-degInG[temp_node];
                povit=temp_node;
            }
        }

        for(temp_i=Cand_R.vnum-1;temp_i>=0;--temp_i){
            temp_node=Cand_R.vlist[temp_i];
            if(temp_j-degInG[temp_node]>temp_value){
                temp_value=temp_j-degInG[temp_node];
                povit=temp_node;
            }
        }

        if(temp_value<=k){
            if(DEBUG) cout<<"not find pivot"<<endl;
            if(Exc_L.empty()&&Exc_R.empty()){
                ++total_num;

                temp_value=0;
                for(temp_i=P_L.vnum-1;temp_i>=0;--temp_i){
                    temp_node=P_L.vlist[temp_i];
                    temp_value+=degInG[temp_node];
                }
                for(temp_i=Cand_L.vnum-1;temp_i>=0;--temp_i){
                    temp_node=Cand_L.vlist[temp_i];
                    temp_value+=degInG[temp_node];
                }
                if(edges<temp_value){
                    // edges=temp_value;
                    // v_l=P_L.vnum+Cand_L.vnum;
                    // v_r=P_R.vnum+Cand_R.vnum;


                    //#ifdef OUTPUT
                    // maxbp->edges=temp_value;
                    // maxbp->left_size=v_l;
                    // maxbp->right_size=v_r;
                    // maxbp->left.clear();
                    // maxbp->right.clear();
                    // for(int ptr_i=P_L.vnum-1;ptr_i>=0;--ptr_i)  maxbp->left.push_back(re_temp_index[P_L.vlist[ptr_i]]);
                    // for(int ptr_i=Cand_L.vnum-1;ptr_i>=0;--ptr_i)  maxbp->left.push_back(re_temp_index[Cand_L.vlist[ptr_i]]);
                    // for(int ptr_i=P_R.vnum-1;ptr_i>=0;--ptr_i)  maxbp->right.push_back(re_temp_index[P_R.vlist[ptr_i]]);
                    // for(int ptr_i=Cand_R.vnum-1;ptr_i>=0;--ptr_i)  maxbp->right.push_back(re_temp_index[Cand_R.vlist[ptr_i]]);
                    //#endif
                    if(DEBUG){
                    cout<<"find a biplex and size is "<<temp_value<<" and now maximum is "<<edges<<endl;
                    cout<<"P_L("<<P_L.vnum<<"):";
                    for(int i=0;i<P_L.vnum;i++)
                        cout<<P_L.vlist[i]<<' ';
                    cout<<endl;
                    cout<<"P_R("<<P_R.vnum<<"):";
                    for(int i=0;i<P_R.vnum;i++)
                        cout<<P_R.vlist[i]<<' ';
                    cout<<endl;
                    cout<<"C_L("<<Cand_L.vnum<<"):";
                    for(int i=0;i<Cand_L.vnum;i++)
                        cout<<Cand_L.vlist[i]<<' ';
                    cout<<endl;
                    cout<<"C_R("<<Cand_R.vnum<<"):";
                    for(int i=0;i<Cand_R.vnum;i++)
                        cout<<Cand_R.vlist[i]<<' ';
                    cout<<endl;}
                    if(top_k.size()<topK){
                        subgraph temp_sub(temp_value,P_L.vnum+Cand_L.vnum,P_R.vnum+Cand_R.vnum);
                        Record(temp_sub.left,temp_sub.right);
                        top_k.push(temp_sub);
                    //edges=edges<temp_j&&edges>0?edges:temp_j;
                    }else if(edges<temp_value){
                        subgraph temp_sub(temp_value,P_L.vnum+Cand_L.vnum,P_R.vnum+Cand_R.vnum);
                        Record(temp_sub.left,temp_sub.right);
                        top_k.push(temp_sub);
                        top_k.pop();
                        edges=top_k.top().edges;
                    }
                }
                return;
            }

            for(temp_i=Exc_R.vnum-1;temp_i>=0;--temp_i){
                temp_node=Exc_R.vlist[temp_i];
                if(degInG[temp_node]==temp_j) return;
            }

            for(temp_i=Exc_L.vnum-1;temp_i>=0;--temp_i){
                temp_node=Exc_L.vlist[temp_i];
                if(degInG[temp_node]==temp_k) return;
            }

            temp_vector.clear();
            for(temp_i=P_L.vnum-1;temp_i>=0;--temp_i){
                temp_node=P_L.vlist[temp_i];
                if(degInG[temp_node]==temp_k-k){
                    temp_vector.push_back(temp_node);
                }
            }
            for(temp_i=Cand_L.vnum-1;temp_i>=0;--temp_i){
                temp_node=Cand_L.vlist[temp_i];
                if(degInG[temp_node]==temp_k-k){
                    temp_vector.push_back(temp_node);
                }
            }
            for(temp_i=Exc_R.vnum-1;temp_i>=0;--temp_i){
                temp_node=Exc_R.vlist[temp_i];
                if(degInG[temp_node]>=temp_j-k){
                    temp_value=1;
                    for(int i:temp_vector){
                        if(!bit_G[temp_node].test(i)){
                            temp_value=0;
                            break;
                        }
                    }
                    if(temp_value) return;
                }
            }


            temp_vector.clear();
            for(temp_i=P_R.vnum-1;temp_i>=0;--temp_i){
                temp_node=P_R.vlist[temp_i];
                if(degInG[temp_node]==temp_j-k){
                    temp_vector.push_back(temp_node);
                }
            }
            for(temp_i=Cand_R.vnum-1;temp_i>=0;--temp_i){
                temp_node=Cand_R.vlist[temp_i];
                if(degInG[temp_node]==temp_j-k){
                    temp_vector.push_back(temp_node);
                }
            }
            for(temp_i=Exc_L.vnum-1;temp_i>=0;--temp_i){
                temp_node=Exc_L.vlist[temp_i];
                if(degInG[temp_node]>=temp_k-k){
                    temp_value=1;
                    for(int i:temp_vector){
                        if(!bit_G[temp_node].test(i)){
                            temp_value=0;
                            break;
                        }
                    }
                    if(temp_value) return;
                }
            }


            temp_value=0;
            for(temp_i=P_L.vnum-1;temp_i>=0;--temp_i){
                temp_node=P_L.vlist[temp_i];
                temp_value+=degInG[temp_node];
            }
            for(temp_i=Cand_L.vnum-1;temp_i>=0;--temp_i){
                temp_node=Cand_L.vlist[temp_i];
                temp_value+=degInG[temp_node];
            }
            if(edges<temp_value){
                // edges=temp_value;
                // v_l=P_L.vnum+Cand_L.vnum;
                // v_r=P_R.vnum+Cand_R.vnum;

                //#ifdef OUTPUT
                // maxbp->edges=temp_value;
                // maxbp->left_size=v_l;
                // maxbp->right_size=v_r;
                // maxbp->left.clear();
                // maxbp->right.clear();
                // for(int ptr_i=P_L.vnum-1;ptr_i>=0;--ptr_i)  maxbp->left.push_back(re_temp_index[P_L.vlist[ptr_i]]);
                // for(int ptr_i=Cand_L.vnum-1;ptr_i>=0;--ptr_i)  maxbp->left.push_back(re_temp_index[Cand_L.vlist[ptr_i]]);
                // for(int ptr_i=P_R.vnum-1;ptr_i>=0;--ptr_i)  maxbp->right.push_back(re_temp_index[P_R.vlist[ptr_i]]);
                // for(int ptr_i=Cand_R.vnum-1;ptr_i>=0;--ptr_i)  maxbp->right.push_back(re_temp_index[Cand_R.vlist[ptr_i]]);
                //#endif
                if(DEBUG){
                cout<<"find a biplex and size is "<<temp_value<<" and now maximum is "<<edges<<endl;
                cout<<"P_L("<<P_L.vnum<<"):";
                for(int i=0;i<P_L.vnum;i++)
                    cout<<P_L.vlist[i]<<' ';
                cout<<endl;
                cout<<"P_R("<<P_R.vnum<<"):";
                for(int i=0;i<P_R.vnum;i++)
                    cout<<P_R.vlist[i]<<' ';
                cout<<endl;
                cout<<"C_L("<<Cand_L.vnum<<"):";
                for(int i=0;i<Cand_L.vnum;i++)
                    cout<<Cand_L.vlist[i]<<' ';
                cout<<endl;
                cout<<"C_R("<<Cand_R.vnum<<"):";
                for(int i=0;i<Cand_R.vnum;i++)
                    cout<<Cand_R.vlist[i]<<' ';
                cout<<endl;}
                if(top_k.size()<topK){
                    subgraph temp_sub(temp_value,P_L.vnum+Cand_L.vnum,P_R.vnum+Cand_R.vnum);
                    Record(temp_sub.left,temp_sub.right);
                    top_k.push(temp_sub);
                    //edges=edges<temp_value&&edges>0?edges:temp_value;
                }else if(edges<temp_value){
                    subgraph temp_sub(temp_value,P_L.vnum+Cand_L.vnum,P_R.vnum+Cand_R.vnum);
                    Record(temp_sub.left,temp_sub.right);
                    top_k.push(temp_sub);
                    top_k.pop();
                    edges=top_k.top().edges;
                }
            }

            ++total_num;
            return;
        }

        if(DEBUG) cout<<"find a pivot in C : "<<povit<<endl;
        RandList *Cand1, *Cand2, *P1, *P2, *Exc1, *Exc2;
        int temp_stand1, temp_stand2;
        if(povit<bi_index){
            Cand1=&Cand_L;
            Cand2=&Cand_R;
            P1=&P_L;
            P2=&P_R;
            Exc1=&Exc_L;
            Exc2=&Exc_R;
            temp_stand1=left;
            temp_stand2=right;
        }else{
            Cand1=&Cand_R;
            Cand2=&Cand_L;
            P1=&P_R;
            P2=&P_L;
            Exc1=&Exc_R;
            Exc2=&Exc_L;
            temp_stand1=right;
            temp_stand2=left;
        }
        Cand1->remove(povit);
        Exc1->add(povit);
        for(temp_i=degree[povit]-1;temp_i>=0;--temp_i){
            temp_node=Graph[povit][temp_i];
            --degInG[temp_node];
        }
        
        if(CanBranch())
            {BKmb_Rec();}

        Exc1->remove(povit);
        P1->add(povit);
        for(temp_i=degree[povit]-1;temp_i>=0;--temp_i){
            temp_node=Graph[povit][temp_i];
            degInG[temp_node]++;
            degInP[temp_node]++;
        }

        vector<int> remove_C1, remove_E1, remove_C2, remove_E2;
        remove_C1.reserve(Cand1->vnum);
        remove_C2.reserve(Cand2->vnum);
        remove_E1.reserve(Exc1->vnum);
        remove_E2.reserve(Exc2->vnum);

        if(degInP[povit]==P2->vnum-k){
            for(temp_i=Cand2->vnum-1;temp_i>=0;--temp_i){
                temp_node=Cand2->vlist[temp_i];
                if(!bit_G[povit].test(temp_node)||degInG[temp_node]<temp_stand1-k){
                    Cand2->remove(temp_node);
                    remove_C2.push_back(temp_node);
                    for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                        temp_node2=Graph[temp_node][temp_j];
                        --degInG[temp_node2];
                    }
                }
            }
            for(temp_i=Exc2->vnum-1;temp_i>=0;--temp_i){
                temp_node=Exc2->vlist[temp_i];
                if(!bit_G[povit].test(temp_node)||degInG[temp_node]<temp_stand1-k){
                    Exc2->remove(temp_node);
                    remove_E2.push_back(temp_node);
                }
            }
        }else{
            temp_size1=P1->vnum-k;
            for(temp_i=Cand2->vnum-1;temp_i>=0;--temp_i){
                temp_node=Cand2->vlist[temp_i];
                if(degInP[temp_node]<temp_size1||degInG[temp_node]<temp_stand1-k){
                    Cand2->remove(temp_node);
                    remove_C2.push_back(temp_node);
                    for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                        temp_node2=Graph[temp_node][temp_j];
                        --degInG[temp_node2];
                    }
                }
            }

            for(temp_i=Exc2->vnum-1;temp_i>=0;--temp_i){
                temp_node=Exc2->vlist[temp_i];
                if(degInP[temp_node]<temp_size1||degInG[temp_node]<temp_stand1-k){
                    Exc2->remove(temp_node);
                    remove_E2.push_back(temp_node);
                }
            }
        }

        
        for(temp_i=degree[povit]-1;temp_i>=0;--temp_i){
            temp_node=Graph[povit][temp_i];
            if(Cand2->contains(temp_node)||P2->contains(temp_node))
                G_index[temp_node]=1;
        }
        temp_vector.clear();
        temp_size1=P1->vnum;
        temp_size2=P2->vnum;
        for(temp_i=P2->vnum-1;temp_i>=0;--temp_i){
            temp_node=P2->vlist[temp_i];
            if(degInP[temp_node]==temp_size1-k&&!bit_G[povit].test(temp_node)){
                temp_vector.push_back(temp_node);
            }
        }

        if(!temp_vector.empty()&&Cand2->vnum+P2->vnum>=temp_stand2){
            for(temp_i=Cand1->vnum-1;temp_i>=0&&Cand1->vnum+P1->vnum>=temp_stand1;--temp_i){
                temp_node=Cand1->vlist[temp_i];
                if(degInP[temp_node]==temp_size2) continue;
                if(degInP[temp_node]<temp_size2-k||degInG[temp_node]<temp_stand2-k){
                    Cand1->remove(temp_node);
                    remove_C1.push_back(temp_node);
                    for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                        temp_node2=Graph[temp_node][temp_j];
                        --degInG[temp_node2];
                    }
                }else{
                    for(int i:temp_vector){
                        if(!bit_G[temp_node].test(i)){
                            Cand1->remove(temp_node);
                            remove_C1.push_back(temp_node);
                            for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                                temp_node2=Graph[temp_node][temp_j];
                                --degInG[temp_node2];
                            }
                            break;
                        }
                    }
                    if(Cand1->contains(temp_node)){
                        count=0;
                        for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                            if(G_index[Graph[temp_node][temp_j]]){
                                count++;
                            }
                        }
                        if(count<temp_stand2-2*k){
                            Cand1->remove(temp_node);
                            remove_C1.push_back(temp_node);
                            for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                                temp_node2=Graph[temp_node][temp_j];
                                --degInG[temp_node2];
                            }
                        }
                    }
                }
            }

            for(temp_i=Exc1->vnum-1;temp_i>=0;--temp_i){
                temp_node=Exc1->vlist[temp_i];
                if(degInP[temp_node]==temp_size2) continue;
                if(degInP[temp_node]<temp_size2-k||degInG[temp_node]<temp_stand2-k){
                    Exc1->remove(temp_node);
                    remove_E1.push_back(temp_node);
                }else{
                    for(int i:temp_vector){
                        if(!bit_G[temp_node].test(i)){
                            Exc1->remove(temp_node);
                            remove_E1.push_back(temp_node);
                            break;
                        }
                    }
                }
                
            }
        }
        for(temp_i=degree[povit]-1;temp_i>=0;--temp_i){
            G_index[Graph[povit][temp_i]]=0;
        }
        if(CanBranch()&&P_R.vnum+Cand_R.vnum>=right&&P_L.vnum+Cand_L.vnum>=left)
             {BKmb_Rec();}

        for(int i:remove_C1){
            Cand1->add(i);
            for(temp_i=degree[i]-1;temp_i>=0;--temp_i){
                ++degInG[Graph[i][temp_i]];
            }
        }

        for(int i:remove_C2){
            Cand2->add(i);
            for(temp_i=degree[i]-1;temp_i>=0;--temp_i){
                ++degInG[Graph[i][temp_i]];
            }
        }

        P1->remove(povit);
        Cand1->add(povit);
        for(temp_i=degree[povit]-1;temp_i>=0;--temp_i){
            --degInP[Graph[povit][temp_i]];
        }

        for(int i:remove_E1){
            Exc1->add(i);
        }
        for(int i:remove_E2){
            Exc2->add(i);
        }

        

    }

}

bool BKMB::CanBranch(){
    temp_j=P_L.vnum+Cand_L.vnum;
    temp_k=P_R.vnum+Cand_R.vnum;
    for(temp_i=Exc_L.vnum-1;temp_i>=0;--temp_i){
        temp_node=Exc_L.vlist[temp_i];
        if(degInG[temp_node]==temp_k){
            return false;
        }
    }

    for(temp_i=Exc_R.vnum-1;temp_i>=0;--temp_i){
        temp_node=Exc_R.vlist[temp_i];
        if(degInG[temp_node]==temp_j){
            return false;
        }
    }

    return true;
}