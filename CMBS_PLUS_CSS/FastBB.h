#pragma once
#include<iostream>
#include<vector>
#include<time.h>
#include"RandList.h"
#include"subgraph.h"
#include<queue>
using namespace std;

class BKMB
{
private:
    int **Graph;
    int *degree;
    int graph_size;
    int bi_index;
    int k;
    int left, right;
    int temp_i,temp_j,temp_k, temp_value, temp_node, temp_node2, temp_size1, temp_size2;
    MBitSet *bit_G;
    RandList P_L, P_R, Cand_L, Cand_R;
    int *degInP;
    int *degInG;
    int *G_index; int count;
    vector<int> temp_vector;
    bool CanBranch();
    void BKmb_Rec();
    void BKmb_Rec2();
    void Record1(vector<int>& left, vector<int>& right);
    void Record2(vector<int>& left, vector<int>& right);
public:
    BKMB(int **Graph, int *degree, int graph_size, int bi_index, int k,int left, int right);
    ~BKMB();
    void BKmb();
    void BKmb2();
    int total_num;
    int edges,v_l,v_r;
    int myset;
    subgraph maxbp;
    subgraph bp;
    bool exist;
};

BKMB::BKMB(int **Graph, int *degree, int graph_size, int bi_index, int k, int left, int right)
{
    this->Graph=Graph;
    this->degree=degree;
    this->graph_size=graph_size;
    this->bi_index=bi_index;
    this->k=k;
    this->left=left;
    this->right=right;
    total_num=0;
    myset=graph_size;
    P_L.init(graph_size);
    P_R.init(graph_size);
    Cand_L.init(graph_size);
    Cand_R.init(graph_size);
    temp_vector.reserve(graph_size);
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
    maxbp=subgraph(0,0,0);
    bp=subgraph(0,0,0);
    exist=false;
}
BKMB::~BKMB()
{
    delete [] degInG;
    delete [] degInP;
    delete [] G_index;
}
void BKMB::BKmb(){
    bool flag;
    for(int i=0;i<bi_index;++i){
        if(edges>=left*right) break;
    // each node i in L to select B_i
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

        BKmb_Rec();
        // clear all set
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
    }
}

void BKMB::BKmb2(){
    bool flag;
    for(int i=0;i<bi_index;++i){
        if(edges>=left*right) break;
        // each node i in L to select B_i
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
        // cout<<"go to branch "<<endl;
        BKmb_Rec2();
        // clear all set
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
    }
}
void BKMB::BKmb_Rec(){
    if(edges>=left*right) return;
    //myset is the upper bound of the L side 
    if(P_L.vnum>=myset){
        return;
    }// if |L|> upper_l return
    if(P_L.vnum+Cand_L.vnum<left||P_R.vnum+Cand_R.vnum<right) return;
    //if |S_L+C_L|<low_l or |S_R+C_R| < low bound
    //if |S|>upper bound x,y
    if(P_L.vnum>=left&&P_R.vnum>=right)
    {
        total_num++;
        temp_j=0;
        for(temp_i=P_L.vnum-1;temp_i>=0;--temp_i){
            temp_node=P_L.vlist[temp_i];
            temp_j+=degInP[temp_node];
        }
        //?when now edges > current edges
        if(edges<temp_j){
            subgraph temp_sub(temp_j,P_L.vnum,P_R.vnum);
            maxbp=temp_sub;Record1(maxbp.left,maxbp.right);
            edges=temp_j;
            if(edges>=left*right) return;    
        }
        return;
    }   
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
    int lowbound=(left<right)?(left*right-left*k):(left*right-right*k);
    if(temp_j<lowbound)
    {
        return;
    }
    // if there don't have nodes in Candidate sets
    if(Cand_R.empty()&&Cand_L.empty()){
        if(P_L.vnum==left&&P_R.vnum==right){
            total_num++;
            temp_j=0;
            for(temp_i=P_L.vnum-1;temp_i>=0;--temp_i){
                temp_node=P_L.vlist[temp_i];
                temp_j+=degInP[temp_node];
            }
            //?when now edges > current edges
            if(edges<temp_j){
                subgraph temp_sub(temp_j,P_L.vnum+Cand_L.vnum,P_R.vnum+Cand_R.vnum);
                maxbp=temp_sub;Record2(maxbp.left,maxbp.right);
                edges=temp_j;
                if(edges>=left*right) return;    
            }
        } 
        return;
    }

    // cout<<"judge finished! now find a povit to order the branch"<<endl;
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
    
    // cout<<"now find a povit: "<<povit<<endl;
    if(povit>=0&&temp_value>k){
        // cout<<"can be povit-based branch "<<endl;
        RandList *Cand1, *Cand2, *P1, *P2;
        int temp_stand1, temp_stand2;
        if(povit<bi_index){
            Cand1=&Cand_L;
            Cand2=&Cand_R;
            P1=&P_L;
            P2=&P_R;
            temp_stand1=left;
            temp_stand2=right;
        }else{
            Cand1=&Cand_R;
            Cand2=&Cand_L;
            P1=&P_R;
            P2=&P_L;
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
        
        for(temp_i=degree[rec_a]-1;temp_i>=0;--temp_i){
            --degInG[Graph[rec_a][temp_i]];
        }
        //first branch  s+=NULL,c remove first element,d add first element
        if(CanBranch()&&P2->vnum+Cand2->vnum>=temp_stand2)
            BKmb_Rec(); //if canbranch: there is no node n in D,make branch not maximum
        if(edges>=left*right) return;
        vector<int> remove_C1, remove_E1, remove_C2, remove_E2;
        remove_C1.reserve(Cand1->vnum); 
        remove_C2.reserve(Cand2->vnum);
      
        //while add node which in gamma(povit,cand) to branch don't violate the k-biplex definition
        while(idx<p+1){
            //except the first branch, other branches contain the first node
            P2->add(rec_a);
            if(P2->vnum>temp_stand2) {/*cout<<"warninng!!!!!!!P is over the upper bound of size!"<<endl; */break;}
            for(temp_i=degree[rec_a]-1;temp_i>=0;--temp_i){
                temp_node=Graph[rec_a][temp_i];
                degInG[temp_node]++;
                degInP[temp_node]++;
            }
            //judged by degree, pruning the C and E after adding the rec_a
            if(degInP[rec_a]==P1->vnum-k){ //if added rec_a and gamma(rec_a,S)=k: any nodes disconnect to rec_a are prunned
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
                }//remove the node don't connect to rec_a or deg don't satisfy the constraint,if remove update the degree in G
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
                }//remove the node d(n)<|P2|-k or don't satisfy the constraint
            }
            
            if(Cand1->vnum+P1->vnum<temp_stand1){
                break;
            }//if the branch don't satisfy the constraint, branches after now branche are prunned

            temp_vector.clear();
            temp_size1=P2->vnum;
            temp_size2=P1->vnum;
            for(temp_i=P1->vnum-1;temp_i>=0;--temp_i){
                temp_node=P1->vlist[temp_i];
                if(degInP[temp_node]==temp_size1-k&&!bit_G[rec_a].test(temp_node)){
                    temp_vector.push_back(temp_node);
                }
            }//temp_vector store the node disconnect to rec_a in S but can not accept the disconnect
            for(temp_i=degree[rec_a]-1;temp_i>=0;--temp_i){
                temp_node=Graph[rec_a][temp_i];
                if(Cand1->contains(temp_node)||P1->contains(temp_node))
                    G_index[temp_node]=1;
            }//G-index: N(rec_a) also in (S or C)


            if(!temp_vector.empty()){
                for(temp_i=Cand2->vnum-1;temp_i>=0;--temp_i){
                    temp_node=Cand2->vlist[temp_i];
                    if(degInP[temp_node]==temp_size2) continue; //if node n in Cand2 all connect to P1,it must be added into P1,continue
                    if(degInP[temp_node]<temp_size2-k||degInG[temp_node]<temp_stand1-k){//if d(n,S)<|P1|-k or d(n,S and C)<theta-k,prune
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
                        }//if can not connect to n:n in S and d(n,S)=|P1|-k
                        
                        //theta-2k constarint to rec_a
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
            }

            for(temp_i=degree[rec_a]-1;temp_i>=0;--temp_i){
                temp_node=Graph[rec_a][temp_i];
                G_index[temp_node]=0;
            }
            

            if(Cand2->vnum+P2->vnum<temp_stand2){
                break;
            }//prune if don't satisfy the constraint

            //find a node in cand2 and doing,p means the number of nodes can be add into S in doing/N(pivot)
            for(;idx<doing.size()&&!Cand2->contains(doing[idx]);++p,++idx);
            if(idx==doing.size()) break;
            rec_b=doing[idx];


            //branch
            Cand2->remove(rec_b);
            for(temp_i=degree[rec_b]-1;temp_i>=0;--temp_i){
                temp_node=Graph[rec_b][temp_i];
                --degInG[temp_node];
            }
            if(CanBranch())
                BKmb_Rec();
            if(edges>=left*right) return;
            rec_a=rec_b;
            ++idx;
        }


        //remove the remain nodes in doing
        //if S+C satisfy the theta constraint
        if(Cand1->vnum+P1->vnum>=temp_stand1&&Cand2->vnum+P2->vnum>=temp_stand2){
            
            count=0;
            for(;idx<doing.size();++idx){
                temp_node=doing[idx];
                if(Cand2->contains(temp_node)){
                    count++;
                }
            }//count : nodes number in doing and Cand
            
            //if S+C-count also satisfy the constraint
            if(Cand2->vnum+P2->vnum-count>=temp_stand2){
                for(;idx<doing.size();++idx){
                temp_node=doing[idx];
                if(Cand2->contains(temp_node)){
                    Cand2->remove(temp_node);
                    for(temp_i=degree[temp_node]-1;temp_i>=0;--temp_i){
                        --degInG[Graph[temp_node][temp_i]];
                    }
                }
                }//remove
                if(CanBranch())
                    BKmb_Rec();
                if(edges>=left*right) return;
            }
            
            
            
        }
        //recall
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

    }else{
        // cout<<"can not povit-based branch"<<endl;
        //if can't find a povit in S find a povit in C
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

        //if not has a povit condition
        if(temp_value<=k){
            // cout<<"can not find a povit in Cand "<<endl;
            if(P_L.vnum+Cand_L.vnum<left||P_R.vnum+Cand_R.vnum<right) 
            {
                return;
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
            cout<<"now construct a biplex  U: "<<P_L.vnum+Cand_L.vnum<<" V: "<<P_R.vnum+Cand_R.vnum<<" size: "<<temp_value<<endl;
            if(edges<temp_value){                
                subgraph temp_sub(temp_value,P_L.vnum+Cand_L.vnum,P_R.vnum+Cand_R.vnum);
                maxbp=temp_sub;edges=temp_value;Record2(maxbp.left,maxbp.right);
                if(edges>=left*right) return; 
            }
            ++total_num;
            return;
        }        
        // cout<<"find a povit in Cand : "<<povit<<endl;
        RandList *Cand1, *Cand2, *P1, *P2;
        int temp_stand1, temp_stand2;
        if(povit<bi_index){
            Cand1=&Cand_L;
            Cand2=&Cand_R;
            P1=&P_L;
            P2=&P_R;
            temp_stand1=left;
            temp_stand2=right;
        }else{
            Cand1=&Cand_R;
            Cand2=&Cand_L;
            P1=&P_R;
            P2=&P_L;
            temp_stand1=right;
            temp_stand2=left;
        }

        Cand1->remove(povit);
        for(temp_i=degree[povit]-1;temp_i>=0;--temp_i){
            temp_node=Graph[povit][temp_i];
            --degInG[temp_node];
        }
        
        if(CanBranch())
            BKmb_Rec();
        if(edges>=left*right) return;

        P1->add(povit);
        for(temp_i=degree[povit]-1;temp_i>=0;--temp_i){
            temp_node=Graph[povit][temp_i];
            degInG[temp_node]++;
            degInP[temp_node]++;
        }

        
        vector<int> remove_C1, remove_E1, remove_C2, remove_E2;
        remove_C1.reserve(Cand1->vnum);
        remove_C2.reserve(Cand2->vnum);

        
        
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
        }
        for(temp_i=degree[povit]-1;temp_i>=0;--temp_i){
            G_index[Graph[povit][temp_i]]=0;
        }
       
        if(CanBranch()&&P_R.vnum+Cand_R.vnum>=right&&P_L.vnum+Cand_L.vnum>=left)
            BKmb_Rec();
        if(edges>=left*right) return;
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
    }
}

void BKMB::BKmb_Rec2(){
    // cout<<"now P_L: ";
    // for(int i=0;i<P_L.vnum;i++)
    // {
    //     cout<<P_L.vlist[i]<<" ";
    // }
    // cout<<endl;
    // cout<<"now P_R: ";
    // for(int i=0;i<P_R.vnum;i++)
    // {
    //     cout<<P_R.vlist[i]<<" ";
    // }
    // cout<<endl;
    // cout<<"now Cand_L: ";
    // for(int i=0;i<Cand_L.vnum;i++)
    // {
    //     cout<<Cand_L.vlist[i]<<" ";
    // }
    // cout<<endl;
    // cout<<"now Cand_R: ";
    // for(int i=0;i<Cand_R.vnum;i++)
    // {
    //     cout<<Cand_R.vlist[i]<<" ";
    // }
    // cout<<endl;
    // cout<<"constraint judge"<<endl;
    if(edges>=left*right) return; 
    if(P_L.vnum>=myset){
        return;
    }
    if(P_L.vnum+Cand_L.vnum<left||P_R.vnum+Cand_R.vnum<right) return;
    if(P_L.vnum>=left&&P_R.vnum>=right)
    {
        cout<<" over the size"<<endl;
        total_num++;
        temp_j=0;
        for(temp_i=P_L.vnum-1;temp_i>=0;--temp_i){
            temp_node=P_L.vlist[temp_i];
            temp_j+=degInP[temp_node];
        }
        subgraph temp_sub(temp_j,P_L.vnum,P_R.vnum);
        if(bp.edges<temp_j&&!exist)
        {
            bp=temp_sub;
            exist=true;
        }
        //?when now edges > current edges
        if(edges<temp_j){
            maxbp=temp_sub;Record1(maxbp.left,maxbp.right);
            edges=temp_j;
            if(edges>=left*right) return;    
        }
        return;
    } 
    if(exist) return;   
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
    if(temp_j<edges&&exist)
    {
        return;
    }
    int lowbound=(left<right)?(left*right-left*k):(left*right-right*k);
    if(temp_j<lowbound)
    {
        return;
    }
    
    // if there don't have nodes in Candidate sets
    if(Cand_R.empty()&&Cand_L.empty()){
        if(P_L.vnum==left&&P_R.vnum==right){
            total_num++;
            temp_j=0;
            for(temp_i=P_L.vnum-1;temp_i>=0;--temp_i){
                temp_node=P_L.vlist[temp_i];
                temp_j+=degInP[temp_node];
            }
            //?when now edges > current edges
            subgraph temp_sub(temp_j,P_L.vnum+Cand_L.vnum,P_R.vnum+Cand_R.vnum);
            if(bp.edges<temp_j&&!exist)
            {
                bp=temp_sub;
                exist=true;
            }
            if(edges<temp_j){
                
                maxbp=temp_sub;Record2(maxbp.left,maxbp.right);
                edges=temp_j;
                if(edges>=left*right) return;    
            }
        } 
        return;
    }

    // cout<<"judge finished! now find a povit to order the branch"<<endl;
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
    
    // cout<<"now find a povit: "<<povit<<endl;
    if(povit>=0&&temp_value>k){
        // cout<<"can be povit-based branch "<<endl;
        RandList *Cand1, *Cand2, *P1, *P2;
        int temp_stand1, temp_stand2;
        if(povit<bi_index){
            Cand1=&Cand_L;
            Cand2=&Cand_R;
            P1=&P_L;
            P2=&P_R;
            temp_stand1=left;
            temp_stand2=right;
        }else{
            Cand1=&Cand_R;
            Cand2=&Cand_L;
            P1=&P_R;
            P2=&P_L;
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
        
        for(temp_i=degree[rec_a]-1;temp_i>=0;--temp_i){
            --degInG[Graph[rec_a][temp_i]];
        }
        //first branch  s+=NULL,c remove first element,d add first element
        if(CanBranch()&&P2->vnum+Cand2->vnum>=temp_stand2)
            BKmb_Rec2(); //if canbranch: there is no node n in D,make branch not maximum
        if(edges>=left*right) return;
        vector<int> remove_C1, remove_E1, remove_C2, remove_E2;
        remove_C1.reserve(Cand1->vnum); 
        remove_C2.reserve(Cand2->vnum);
      
        //while add node which in gamma(povit,cand) to branch don't violate the k-biplex definition
        while(idx<p+1){
            //except the first branch, other branches contain the first node
            P2->add(rec_a);
            if(P2->vnum>temp_stand2) {/*cout<<"warninng!!!!!!!P is over the upper bound of size!"<<endl; */break;}
            for(temp_i=degree[rec_a]-1;temp_i>=0;--temp_i){
                temp_node=Graph[rec_a][temp_i];
                degInG[temp_node]++;
                degInP[temp_node]++;
            }
            //judged by degree, pruning the C and E after adding the rec_a
            if(degInP[rec_a]==P1->vnum-k){ //if added rec_a and gamma(rec_a,S)=k: any nodes disconnect to rec_a are prunned
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
                }//remove the node don't connect to rec_a or deg don't satisfy the constraint,if remove update the degree in G
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
                }//remove the node d(n)<|P2|-k or don't satisfy the constraint
            }
            
            if(Cand1->vnum+P1->vnum<temp_stand1){
                break;
            }//if the branch don't satisfy the constraint, branches after now branche are prunned

            temp_vector.clear();
            temp_size1=P2->vnum;
            temp_size2=P1->vnum;
            for(temp_i=P1->vnum-1;temp_i>=0;--temp_i){
                temp_node=P1->vlist[temp_i];
                if(degInP[temp_node]==temp_size1-k&&!bit_G[rec_a].test(temp_node)){
                    temp_vector.push_back(temp_node);
                }
            }//temp_vector store the node disconnect to rec_a in S but can not accept the disconnect
            for(temp_i=degree[rec_a]-1;temp_i>=0;--temp_i){
                temp_node=Graph[rec_a][temp_i];
                if(Cand1->contains(temp_node)||P1->contains(temp_node))
                    G_index[temp_node]=1;
            }//G-index: N(rec_a) also in (S or C)


            if(!temp_vector.empty()){
                for(temp_i=Cand2->vnum-1;temp_i>=0;--temp_i){
                    temp_node=Cand2->vlist[temp_i];
                    if(degInP[temp_node]==temp_size2) continue; //if node n in Cand2 all connect to P1,it must be added into P1,continue
                    if(degInP[temp_node]<temp_size2-k||degInG[temp_node]<temp_stand1-k){//if d(n,S)<|P1|-k or d(n,S and C)<theta-k,prune
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
                        }//if can not connect to n:n in S and d(n,S)=|P1|-k
                        
                        //theta-2k constarint to rec_a
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
            }

            for(temp_i=degree[rec_a]-1;temp_i>=0;--temp_i){
                temp_node=Graph[rec_a][temp_i];
                G_index[temp_node]=0;
            }
            

            if(Cand2->vnum+P2->vnum<temp_stand2){
                break;
            }//prune if don't satisfy the constraint

            //find a node in cand2 and doing,p means the number of nodes can be add into S in doing/N(pivot)
            for(;idx<doing.size()&&!Cand2->contains(doing[idx]);++p,++idx);
            if(idx==doing.size()) break;
            rec_b=doing[idx];


            //branch
            Cand2->remove(rec_b);
            for(temp_i=degree[rec_b]-1;temp_i>=0;--temp_i){
                temp_node=Graph[rec_b][temp_i];
                --degInG[temp_node];
            }
            if(CanBranch())
                BKmb_Rec2();
            if(edges>=left*right) return;
            rec_a=rec_b;
            ++idx;
        }


        //remove the remain nodes in doing
        //if S+C satisfy the theta constraint
        if(Cand1->vnum+P1->vnum>=temp_stand1&&Cand2->vnum+P2->vnum>=temp_stand2){
            
            count=0;
            for(;idx<doing.size();++idx){
                temp_node=doing[idx];
                if(Cand2->contains(temp_node)){
                    count++;
                }
            }//count : nodes number in doing and Cand
            
            //if S+C-count also satisfy the constraint
            if(Cand2->vnum+P2->vnum-count>=temp_stand2){
                for(;idx<doing.size();++idx){
                temp_node=doing[idx];
                if(Cand2->contains(temp_node)){
                    Cand2->remove(temp_node);
                    for(temp_i=degree[temp_node]-1;temp_i>=0;--temp_i){
                        --degInG[Graph[temp_node][temp_i]];
                    }
                }
                }//remove
                if(CanBranch())
                    BKmb_Rec2();
                if(edges>=left*right) return;
            }
            
            
            
        }
        //recall
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

    }else{
        // cout<<"can not povit-based branch"<<endl;
        //if can't find a povit in S find a povit in C
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

        //if not has a povit condition
        if(temp_value<=k){
            // cout<<"can not find a povit in Cand "<<endl;
            if(P_L.vnum+Cand_L.vnum<left||P_R.vnum+Cand_R.vnum<right) 
            {
                return;
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
            // cout<<"now construct a biplex  U: "<<P_L.vnum+Cand_L.vnum<<" V: "<<P_R.vnum+Cand_R.vnum<<" size: "<<temp_value<<endl;
            subgraph temp_sub(temp_value,P_L.vnum+Cand_L.vnum,P_R.vnum+Cand_R.vnum);
            if(bp.edges<temp_j&&!exist)
            {
                bp=temp_sub;
                exist=true;
            }
            if(edges<temp_value){                
                maxbp=temp_sub;Record2(maxbp.left,maxbp.right);
                edges=temp_value;
                if(edges>=left*right) return; 
            }
            ++total_num;
            return;
        }        
        // cout<<"find a povit in Cand : "<<povit<<endl;
        RandList *Cand1, *Cand2, *P1, *P2;
        int temp_stand1, temp_stand2;
        if(povit<bi_index){
            Cand1=&Cand_L;
            Cand2=&Cand_R;
            P1=&P_L;
            P2=&P_R;
            temp_stand1=left;
            temp_stand2=right;
        }else{
            Cand1=&Cand_R;
            Cand2=&Cand_L;
            P1=&P_R;
            P2=&P_L;
            temp_stand1=right;
            temp_stand2=left;
        }

        Cand1->remove(povit);
        for(temp_i=degree[povit]-1;temp_i>=0;--temp_i){
            temp_node=Graph[povit][temp_i];
            --degInG[temp_node];
        }
        
        if(CanBranch())
            BKmb_Rec2();
        if(edges>=left*right) return;

        P1->add(povit);
        for(temp_i=degree[povit]-1;temp_i>=0;--temp_i){
            temp_node=Graph[povit][temp_i];
            degInG[temp_node]++;
            degInP[temp_node]++;
        }

        
        vector<int> remove_C1, remove_E1, remove_C2, remove_E2;
        remove_C1.reserve(Cand1->vnum);
        remove_C2.reserve(Cand2->vnum);

        
        
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
        }
        for(temp_i=degree[povit]-1;temp_i>=0;--temp_i){
            G_index[Graph[povit][temp_i]]=0;
        }
       
        if(CanBranch()&&P_R.vnum+Cand_R.vnum>=right&&P_L.vnum+Cand_L.vnum>=left)
            BKmb_Rec2();
        if(edges>=left*right) return;
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
    }
}

bool BKMB::CanBranch(){
    if(P_L.vnum>left||P_R.vnum>right) return false;
    return true;
}
void BKMB::Record1(vector<int>& left, vector<int>& right){
    for(int i=P_L.vnum-1;i>=0;--i)
        left.push_back(P_L.vlist[i]);
    for(int i=P_R.vnum-1;i>=0;--i)
        right.push_back(P_R.vlist[i]);
}
void BKMB::Record2(vector<int>& left, vector<int>& right){
    for(int i=P_L.vnum-1;i>=0;--i)
        left.push_back(P_L.vlist[i]);
    for(int i=Cand_L.vnum-1;i>=0;--i)
        left.push_back(Cand_L.vlist[i]);
    for(int i=P_R.vnum-1;i>=0;--i)
        right.push_back(P_R.vlist[i]);
    for(int i=Cand_R.vnum-1;i>=0;--i)
        right.push_back(Cand_R.vlist[i]);
}