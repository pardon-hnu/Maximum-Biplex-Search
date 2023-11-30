#include "biplex.h"

biplexminer::biplexminer(int **Graph, int Graph_size, int Bipartite_index,int *degree, int S, int k){
    this->Graph=Graph;
    this->Bipartite_index=Bipartite_index;
    this->S=S;
    this->k=k;
    this->Graph_size=Graph_size;
    this->degree=degree;


    G_State=new int[Graph_size];
    G_index=new int[Graph_size];
    G_label=new int[Graph_size];
    G_temp=new int[Graph_size];
    G_prune=new int[Graph_size];


    res_count=0;
    res_size=0;
    
    for(int i=0;i<Graph_size;++i){
        G_index[i]=0;
        G_State[i]=0;
        G_label[i]=0;
        G_temp[i]=0;
        G_prune[i]=0;
    }

    for(int i=Bipartite_index;i<Graph_size;++i){
        G_label[i]=1;
    }

}

biplexminer::~biplexminer(){
    delete []G_State;
    delete []G_index;
    delete []G_label;
    delete []G_temp;
    delete []G_prune;
}

void biplexminer::Miner(){
    vector<int> X;
    vector<int> N_X;
    vector<int> cand_q;
    vector<int> Cand_exts;
    for(int i=0;i<Bipartite_index;++i){
        Cand_exts.push_back(i);
    }
    for (int i = Bipartite_index; i < Graph_size; i++){
        N_X.push_back(i);
    }
    GenerateT=0;
    EnumerateT=0;
    time_t e1=clock();
    Enumerate(X,N_X,Cand_exts,cand_q);
    time_t e2=clock();
    EnumerateT=(double)(e2-e1)/CLOCKS_PER_SEC-GenerateT;
    //cout<<"generate cost: "<<GenerateT<<"s."<<endl;
    //cout<<"enumerate cost: "<<EnumerateT<<"s. "<<endl;

    //cout<<"Now print the enumerate X'set size:"<<Xs.size()<<endl;
    // for(int i=0;i<Xs.size();i++)
    // {
    //     for(int j=0;j<Xs[i].size();j++)
    //     {
    //         //cout<<Xs[i][j]<<" ";
    //     }
    //     //cout<<endl;
    // }
}

void biplexminer::Enumerate(vector<int> &X, vector<int> &N_X, vector<int> &Cand_p, vector<int> &Cand_q){
    Xs.push_back(X);
    int temp_node;
    int temp_size;
    vector<int> new_N_X;
    vector<int> new_X;
    vector<int> new_Cand_p;
    vector<int> new_Cand_q;
    vector<int>:: iterator iter1;
    vector<int>:: iterator iter2;

    // //cout<<"the "<<Xs.size()<<" th enumerate is begin"<<endl;
    // //cout<<"print the four argument "<<endl;
    // //cout<<"X: ";
    // for(int i=0;i<X.size();i++)
    // //cout<<X[i]<<' ';
    // //cout<<endl;
    // //cout<<"N_X: ";
    // for(int i=0;i<N_X.size();i++)
    // //cout<<N_X[i]<<' ';
    // //cout<<endl;
    // //cout<<"Cand_p: ";
    // for(int i=0;i<Cand_p.size();i++)
    // //cout<<Cand_p[i]<<' ';
    // //cout<<endl;
    // //cout<<"Cand_q: ";
    // for(int i=0;i<Cand_q.size();i++)
    // //cout<<Cand_q[i]<<' ';
    // //cout<<endl;
    //X is V, and N_X is N(V).Cand_p is the candidate set of V, Cand_q ?,G_state is the degree, G_label is the node in X and N_X,
    // //cout<<"now copy cand_q to new_cand_q"<<endl;
    for(iter1=Cand_q.begin();iter1!=Cand_q.end();++iter1){
        new_Cand_q.push_back(*iter1);   
    }
    // //cout<<"while Cand_p is not empty"<<endl;
    while(!Cand_p.empty()){

        new_N_X.clear();
        new_X.clear();
        new_Cand_p.clear();

        // //cout<<"copy X to new_X"<<endl;
        for(iter1=X.begin();iter1!=X.end();++iter1){
            new_X.push_back(*iter1);
        }

        //temp_node: v from Cand_V
        int temp_node=Cand_p.front();
        new_X.push_back(temp_node);
        //note G_label(v)=1 and erase v from Cand_V
        G_label[temp_node]=1;
        Cand_p.erase(Cand_p.begin());
        // //cout<<"insert "<<temp_node<<" into new_X ,and delete it from cand_p"<<endl;

        if((int)(new_X.size())<= k+1&&(int)(new_X.size())<S){
            //  //cout<<"|new_X|<=k+1 and |new_X|<S"<<endl;
            //update G_State: in G_state each node u connect to v, u's value+=1
            for(int i=0;i<degree[temp_node];++i){
                G_State[Graph[temp_node][i]]+=1;
            }
            // //cout<<"update G_state to: ";
            // for(int i=0;i<Graph_size;i++)
            // //cout<<i<<":"<<G_State[i]<<' ';
            // //cout<<endl;
            //if new_V size == k+1,insert each u which connect to v into new_N(V) and note else u label=0
            //else copy N(V) to new_N(V)
            if(new_X.size()== k+1){
                // //cout<<"new_X size is k+1 so update new_N_X and G_state"<<endl;
                for(int i=Bipartite_index;i<Graph_size;++i){
                    if(G_State[i]>0)
                        new_N_X.push_back(i);
                    else
                        G_label[i]=0;
                }
            }else{
                // //cout<<"new_X size is not k+1 so copy N_X to new_N_X"<<endl;
                new_N_X=N_X;
            }
            // //cout<<"update new_N_X to: ";
            // for(int i=0;i<new_N_X.size();i++)
            // //cout<<i<<":"<<new_N_X[i]<<' ';
            // //cout<<endl;
            // //cout<<"update G_label to: ";
            // for(int i=0;i<Graph_size;i++)
            // //cout<<i<<":"<<G_label[i]<<' ';
            // //cout<<endl;

            //update G_State[v]=degree[v]
            temp_size=degree[temp_node];
            G_State[temp_node]=temp_size;

            // //cout<<"update G_state to: ";
            // for(int i=0;i<Graph_size;i++)
            // //cout<<i<<":"<<G_State[i]<<' ';
            // //cout<<endl;

            // //cout<<"update the cand_p"<<endl;
            for(iter1=Cand_p.begin();iter1!=Cand_p.end();++iter1){

                // //cout<<"now deal the node "<<*iter1<<" in cand_p"<<endl;
                //if d(v')<S-k continue
                temp_size=degree[*iter1];
                if(temp_size<S-k){
                    // //cout<<"do not satisfy the degree constraint S-k,next node"<<endl;
                    continue;
                }

                //possible debug label
                //foreach u' in N(v')
                //if label(u')==1 index(u')==1
                for(int i=0;i<temp_size;++i){

                    if(G_label[Graph[*iter1][i]]==1)
                        G_index[Graph[*iter1][i]]=1;
                }
                // //cout<<"G_index: ";
                // for(int ni=0;ni<Graph_size;ni++)
                // //cout<<ni<<":"<<G_index[ni]<<' ';
                // //cout<<endl;
                //foreach u in N(v)
                //if index(u)==1 temp_count++
                int temp_count=0;
                for(int i=0;i<degree[temp_node];++i){
                    if(G_index[Graph[temp_node][i]]==1){
                        temp_count++;
                    }
                }
                // //cout<<"temp_count  :  "<<temp_count<<endl;
                //if temp_count<S-2*k 
                //foreach u' in N(v'),reset the index, continue;
                if(temp_count<S-2*k){
                    for(int i=0;i<temp_size;++i){
                        G_index[Graph[*iter1][i]]=0;
                    }
                    continue;
                }
                // //cout<<"insert "<<*iter1<<" into new_can_p"<<endl;
                //if temp_count>=S-2*k,reset the index, then v'--> new_cand_V
                for(int i=0;i<temp_size;++i){
                    G_index[Graph[*iter1][i]]=0;
                }
                new_Cand_p.push_back(*iter1);
            }
            // //cout<<"update new Cand p to: ";
            // for(int ni=0;ni<new_Cand_p.size();ni++)
            // {
            //     //cout<<new_Cand_p[ni]<<' ';
            // }
            // //cout<<endl;
            //check the probability of extension
            if(degree[temp_node]>=S-k){
                bool Extension=true;
                temp_size=degree[temp_node];
                for(int i=0;i<temp_size;++i){
                    //if(G_label[])
                    G_index[Graph[temp_node][i]]=1;
                }
                for(iter1=X.begin();iter1!=X.end();++iter1){
                    int temp_count=0;
                    for(int i=0;i<degree[*iter1];++i){
                        if(G_index[Graph[*iter1][i]]==1){
                            temp_count++;
                        }
                    }
                    if(temp_count<S-2*k){
                        Extension=false;
                        for(int i=0;i<temp_size;++i){
                            G_index[Graph[temp_node][i]]=0;
                        }
                        break;
                    }
                }
                for(int i=0;i<temp_size;++i){
                    G_index[Graph[temp_node][i]]=0;
                }


                if(Extension){
                    Enumerate(new_X,new_N_X,new_Cand_p,new_Cand_q);
                }
            }

            new_Cand_q.push_back(temp_node);
            G_State[temp_node]=0;
            G_label[temp_node]=0;
            for(int i=0;i<degree[temp_node];++i){
                G_State[Graph[temp_node][i]]-=1;
            }
            if(new_X.size()==k+1){
                for(int i=Bipartite_index;i<Graph_size;++i){
                    G_label[i]=1;
                }
            }
        }
        else
        {
            //update degree of Y 
            // //cout<<"satisfy the size constraint"<<endl;
            // //cout<<"update G_state and G_label 's value about v"<<endl;
            bool cost=false;
            int temp_count=0;
            for(int i=0;i<degree[temp_node];++i){
                if(G_label[Graph[temp_node][i]]==1){
                    temp_count++;
                    G_State[Graph[temp_node][i]]+=1;
                }
            }
            //update degree of X
            // //cout<<"now remove some node from N_X and copy N_X to new_N_X judged by degree constraint"<<endl;
            G_State[temp_node]=temp_count;
            for(iter1=N_X.begin();iter1!=N_X.end();++iter1){
                if(G_State[*iter1]>=(int)(new_X.size())-k){
                    new_N_X.push_back(*iter1);
                }else{
                    G_label[*iter1]=0;
                    if(G_State[*iter1]!=0){
                        cost=true;
                        for(int i=0;i<degree[*iter1];++i){
                            if(G_label[Graph[*iter1][i]]==1){
                                G_State[Graph[*iter1][i]]-=1;
                            }
                        }
                    }
                }
            }

            bool Extension=true;
            //pruning 1 : check probability of extension"
            // //cout<<"pruning 1: check each two node in new_X is or not satisfy the S-2k degree constraint"<<endl;
            if(new_N_X.size()<0){
                Extension=false; 
            }
            //cost means N_X has some node didn't satisfy the degree constraint
            if(Extension&&cost){
                bool prune=true;
                for(iter1=new_X.begin();iter1!=new_X.end();++iter1){
                    iter2=iter1;iter2++;
                    temp_size=degree[*iter1];
                    for(int i=0;i<temp_size;++i){
                        if(G_label[Graph[*iter1][i]]==1){
                            G_index[Graph[*iter1][i]]=1;
                        }
                    }
                    for(;iter2!=new_X.end();++iter2){
                        int temp_count=0;
                        for(int i=0;i<degree[*iter2];++i){
                            if(G_index[Graph[*iter2][i]]==1){
                                temp_count++;
                            }
                        }
                        if(temp_count<S-2*k){
                            //prune=false;
                            break;
                        }
                    }

                    for(int i=0;i<temp_size;++i){
                        if(G_label[Graph[*iter1][i]]==1){
                            G_index[Graph[*iter1][i]]=0;
                        }
                    }
                    if(!prune){
                        //Extension=false;
                        break;
                    }
                }

            }

            ////cout<<"prunning1 : "<<Extension<<endl;
            //pruning 2
            // //cout<<"pruning 2: base on the lemma 2 from the yu 2020 TKDE"<<endl;
            if(Extension){
                for(iter1=Cand_q.begin();iter1!=Cand_q.end();++iter1){
                    if(degree[*iter1]>=new_N_X.size()){
                        int temp_count=0;
                        for(int i=0;i<degree[*iter1];++i){
                            if(G_label[Graph[*iter1][i]]==1){
                                ++temp_count;
                            }
                        }
                        if(temp_count>=new_N_X.size()){
                            Extension=false;
                            break;
                        }
                    }
                }
            }
            // //cout<<"prunning2 : "<<Extension<<endl;
            //pruning 3
            vector<int> Component;
            if(Extension){

                for(iter1=Cand_p.begin();iter1!=Cand_p.end();++iter1){
                    if(degree[*iter1]>=S-k){
                        int temp_count=0;
                        for(int i=0;i<degree[*iter1];++i){
                            if(G_label[Graph[*iter1][i]]==1){
                                ++temp_count;
                            }
                        }
                        if(temp_count>=new_N_X.size()){
                            Component.push_back(*iter1);
                            //new_Cand_p.erase(iter1);
                            new_X.push_back(*iter1);
                            G_label[*iter1]=1;
                            G_State[*iter1]=new_N_X.size();
                        }else if(temp_count>=S-k){
                            new_Cand_p.push_back(*iter1);
                        }
                    }
                }
                if(Component.size()!=0){
                    for(iter1=new_N_X.begin();iter1!=new_N_X.end();++iter1){
                        G_State[*iter1]+=(Component.size());
                    }
                }
            }
            ////cout<<"prunning3 : "<<Extension<<endl;
            if(Extension){

                for(iter1=new_Cand_p.begin();iter1!=new_Cand_p.end();++iter1){
                    temp_size=degree[*iter1];
                    int temp_debug=*iter1;
                    for(int i=0;i<temp_size;++i){
                        if(G_label[Graph[*iter1][i]]==1){
                            G_index[Graph[*iter1][i]]=1;
                        }
                    }

                    for(iter2=new_X.begin();iter2!=new_X.end();++iter2){
                        int temp_count=0;
                        for(int i=0;i<degree[*iter2];++i){
                            if(G_index[Graph[*iter2][i]]==1){
                                ++temp_count;
                            }
                        }
                        if(temp_count<S-2*k){
                           iter1= new_Cand_p.erase(iter1);
                            break;
                        }
                    }

                    for(int i=0;i<temp_size;++i){
                        if(G_label[Graph[temp_debug][i]]==1){
                            G_index[Graph[temp_debug][i]]=0;
                        }
                    }

                    if(iter1==new_Cand_p.end())
                        break;
                }

                if(new_X.size()>=S){
                    Generate(new_X,new_N_X,new_Cand_p,new_Cand_q);

                    for(iter1=new_N_X.begin();iter1!=new_N_X.end();++iter1){
                        if(G_prune[*iter1]==0){
                            G_label[*iter1]=-1;

                            for(int i=0;i<degree[*iter1];++i){
                                if(G_label[Graph[*iter1][i]]==1){
                                    G_State[Graph[*iter1][i]]-=1;
                                }
                            }

                            iter1= new_N_X.erase(iter1);
                        }else{
                            G_prune[*iter1]=0;
                        }

                        if(iter1==new_N_X.end()) break;
                    }

                    for(iter1=new_Cand_p.begin();iter1!=new_Cand_p.end();++iter1){
                        if(G_prune[*iter1]==0){
                        }else{
                            G_prune[*iter1]=0;
                        }
                        ////cout<<"######"<<*iter1<<endl;
                        if(iter1==new_Cand_p.end()) break;
                    }
                }

                if(new_X.size()+new_Cand_p.size()>=S&&new_N_X.size()>=S){
                    Enumerate(new_X,new_N_X,new_Cand_p,new_Cand_q);
                }

            }



            new_Cand_q.push_back(temp_node);
            G_State[temp_node]=0;
            G_label[temp_node]=0;
            if(Component.size()!=0){
                for(iter1=Component.begin();iter1!=Component.end();++iter1){
                    G_State[*iter1]=0;
                    G_label[*iter1]=0;
                }
            }

            for(iter1=N_X.begin();iter1!=N_X.end();++iter1){
                if(G_label[*iter1]==1){
                    G_State[*iter1]-=(Component.size());
                }
                if(G_label[*iter1]==-1){
                    G_State[*iter1]-=(Component.size());
                    G_label[*iter1]=0;
                }
                if(G_label[*iter1]==0){
                    G_label[*iter1]=1;
                    if(G_State[*iter1]!=0){
                        for(int i=0;i<degree[*iter1];++i){
                            if(G_label[Graph[*iter1][i]]==1){
                                G_State[Graph[*iter1][i]]+=1;
                            }
                        }
                    }
                }
            }
            for(int i=0;i<degree[temp_node];++i){
                if(G_label[Graph[temp_node][i]]==1)
                    G_State[Graph[temp_node][i]]-=1;
            }

        }
    }

}

void biplexminer::Generate(vector<int> &X, vector<int> &N_X, vector<int> &Cand_p, vector<int> &Cand_q){
    time_t g1=clock();
    vector<int> invalid_X;
    vector<int> constant_Y;
    vector<int> var_Y;
    vector<int> ext_q;
    vector<int> cur_Y;
    vector<int> candY_p;
    vector<int>::iterator iter1;
    vector<int>::iterator iter2;
    int temp_num;

    //find invalid X
    temp_num=(int)(N_X.size())-k;
    for(iter1=X.begin();iter1!=X.end();++iter1){
        if(G_State[*iter1]<temp_num){
            invalid_X.push_back(*iter1);
        }
    }

    ////cout<<"***"<<invalid_X.size()<<endl;
    //find constant_Y
    if(invalid_X.size()>0){
        iter1=invalid_X.begin();
        temp_num=degree[*iter1];
        for(int i=0;i<temp_num;++i){
            if(G_label[Graph[*iter1][i]]==1){
                G_index[Graph[*iter1][i]]=1;
            }
        }
        for(iter1++;iter1!=invalid_X.end();++iter1){
            for(int i=0;i<degree[*iter1];++i){
                if(G_index[Graph[*iter1][i]]>=1){
                    G_index[Graph[*iter1][i]]+=1;
                }
            }
        }
        for(iter1=N_X.begin();iter1!=N_X.end();++iter1){
            if(G_index[*iter1]==invalid_X.size()){
                constant_Y.push_back(*iter1);
            }else{
                var_Y.push_back(*iter1);
            }
            G_index[*iter1]=0;
        }
        iter1=invalid_X.begin();
        for(int i=0;i<temp_num;++i){
            // if(G_index[Graph[*iter1][i]]==invalid_X.size()){
            //     constant_Y.push_back(Graph[*iter1][i]);
            // }
            G_index[Graph[*iter1][i]]=0;
        }


        //degree induced by constant_Y
        for(iter1=constant_Y.begin();iter1!=constant_Y.end();++iter1){
            for(int i=0;i<degree[*iter1];++i){
                if(G_label[Graph[*iter1][i]]==1){
                    G_temp[Graph[*iter1][i]]+=1;
                }
            }
            G_prune[*iter1]=1;
        }
    }else{
        for(iter1=N_X.begin();iter1!=N_X.end();++iter1){
            var_Y.push_back(*iter1);
        }
    }




    //obtain ext_q
    for(iter1=Cand_q.begin();iter1!=Cand_q.end();++iter1){
        int temp_count=0;
        for(int i=0;i<degree[*iter1];++i){
            if(G_label[Graph[*iter1][i]]==1){
                temp_count++;
            }
        }
        if(temp_count>=S-k){
            ext_q.push_back(*iter1);
        }
    }
    ////cout<<"%%%"<<constant_Y.size()<<" "<<var_Y.size()<<endl;
    if(invalid_X.size()==0&&(constant_Y.size()+var_Y.size()>=S)){
        bool res=CheckMaximality(X,constant_Y,var_Y,candY_p,Cand_p,ext_q);
        if(res){
            res_count++;

	       if(res_size<GetSize(X,var_Y,constant_Y))
           {
                res_size = GetSize(X,var_Y,constant_Y);
                Lnode=X;
                Rnode=var_Y;
                for(int ri=0;ri<constant_Y.size();ri++)
                {
                    int curele=constant_Y[ri];
                    Rnode.push_back(curele);
                }
           }
            //Output(X,var_Y,constant_Y);
        }
    }else if(invalid_X.size()==1&&k<=4){
         if((int)constant_Y.size()>=S-k&&(constant_Y.size()+var_Y.size())>=S){
             Special(X,var_Y,constant_Y,cur_Y,candY_p,Cand_p,ext_q);
         }
    }else{
        List(X,var_Y,constant_Y,cur_Y,candY_p,Cand_p,ext_q);
    }
    //recovery
    for(iter1=X.begin();iter1!=X.end();++iter1){
        G_temp[*iter1]=0;
    }
    time_t g2=clock();
    GenerateT+=(double)(g2-g1)/CLOCKS_PER_SEC;
}


void biplexminer::List(vector<int> &X, vector<int> var_Y, vector<int> &con_Y, vector<int> &cur_Y, vector<int> &candY_q, vector<int> &ext_p,
vector<int> &ext_q){
    vector<int>::iterator iter1;
    vector<int>::iterator iter2;
    int temp_node;
    int temp_num;

    vector<int> new_cur_Y;
    vector<int> new_var_Y;
    vector<int> new_CandY_q;

    bool check=true;

    for(iter1=candY_q.begin();iter1!=candY_q.end();++iter1){
        new_CandY_q.push_back(*iter1);
    }

    while(!var_Y.empty()){
        new_cur_Y.clear();
        new_var_Y.clear();


        temp_node=var_Y.front();
        var_Y.erase(var_Y.begin());
        ////cout<<temp_node<<" "<<cur_Y.size()<<endl;

        bool extension=true;
        //update
        for(int i=0;i<degree[temp_node];++i){
            if(G_label[Graph[temp_node][i]]==1){
                G_temp[Graph[temp_node][i]]+=1;
            }
        }
        //check
        temp_num=con_Y.size()+cur_Y.size()+1;
        for(iter1=X.begin();iter1!=X.end();++iter1){
            if(G_temp[*iter1]<temp_num-k){
                extension=false;
                break;
            }
        }

        if(extension){
            check=false;
            for(iter1=cur_Y.begin();iter1!=cur_Y.end();iter1++){
                new_cur_Y.push_back(*iter1);
            }
            //temp_node=var_Y.front();
            new_cur_Y.push_back(temp_node);
            //.pop_front();

            for(iter1=var_Y.begin();iter1!=var_Y.end();iter1++){
                new_var_Y.push_back(*iter1);
            }
            List(X,new_var_Y,con_Y,new_cur_Y,new_CandY_q,ext_p,ext_q);
        }

        //update
        new_CandY_q.push_back(temp_node);
        //recovery
        for(int i=0;i<degree[temp_node];++i){
            if(G_label[Graph[temp_node][i]]==1){
                G_temp[Graph[temp_node][i]]-=1;
            }
        }


    }
    if(var_Y.empty()&&check&&(cur_Y.size()+con_Y.size())>=S){
        bool res=CheckMaximality(X,con_Y,cur_Y,candY_q,ext_p,ext_q);
        if(res){
            res_count++;
	       if(res_size<GetSize(X,cur_Y,con_Y))
           {
                res_size = GetSize(X,cur_Y,con_Y);
                Lnode=X;
                Rnode=cur_Y;
                for(int ri=0;ri<con_Y.size();ri++)
                {
                    int curele=con_Y[ri];
                    Rnode.push_back(curele);
                }
           }
	           
            //Output(X,cur_Y,con_Y);

        }
    }

}

bool biplexminer::CheckMaximality(vector<int> &X, vector<int> &con_Y, vector<int> &cur_Y,
    vector<int> &candY_q, vector<int> &ext_p, vector<int> &ext_q){

    vector<int>::iterator iter1;
    vector<int>::iterator iter2;
    int temp_size;
    int temp_num;

    bool isMax=true;
    temp_num=con_Y.size()+cur_Y.size()+1;

    //check Y
    for(iter1=candY_q.begin();iter1!=candY_q.end();++iter1){
        //update
        bool iteration=true;
        temp_size=degree[*iter1];
        for(int i=0;i<temp_size;++i){
            if(G_label[Graph[*iter1][i]]==1){
                G_temp[Graph[*iter1][i]]+=1;
            }
        }
        //check
        for(iter2=X.begin();iter2!=X.end();++iter2){
            if(G_temp[*iter2]<temp_num-k){
                iteration=false;
                break;
            }
        }
        //recovery
        for(int i=0;i<temp_size;++i){
            if(G_label[Graph[*iter1][i]]==1){
                G_temp[Graph[*iter1][i]]-=1;
            }
        }
        if(iteration){
            return false;
        }
    }
    for(iter1=cur_Y.begin();iter1!=cur_Y.end();++iter1){
            G_index[*iter1]=1;
    }
    for(iter2=con_Y.begin();iter2!=con_Y.end();++iter2){
        G_index[*iter2]=1;
    }

    //optimization
    for(iter1=ext_q.begin();iter1!=ext_q.end();++iter1){
        //mark
        int temp_countt=0;
        bool iteration=true;



        for(int i=0;i<degree[*iter1];++i){
            if(G_index[Graph[*iter1][i]]==1){
                if(isMax){
                    G_index[Graph[*iter1][i]]=2;
                }
                temp_countt++;
            }
        }
        if(temp_countt==temp_num-1){
            ////cout<<*iter1<<" "<<temp_countt<<endl;
            for(iter2=con_Y.begin();iter2!=con_Y.end();++iter2){
                G_index[*iter2]=0;
            }
            for(iter2=cur_Y.begin();iter2!=cur_Y.end();++iter2){
                G_index[*iter2]=0;
            }
            return false;
        }else if(isMax&&temp_countt>=temp_num-1-k){
            for(iter2=con_Y.begin();iter2!=con_Y.end();++iter2){
                if(G_index[*iter2]==1&&G_State[*iter2]<(int)(X.size())+1-k){
                    iteration=false;
                    break;
                }
            }
            for(iter2=cur_Y.begin();iter2!=cur_Y.end();++iter2){
                if(G_index[*iter2]==1&&G_State[*iter2]<(int)(X.size())+1-k){
                    iteration=false;
                    break;
                }
            }
        }else{
            iteration=false;
        }
        if(isMax){
            for(iter2=cur_Y.begin();iter2!=cur_Y.end();++iter2){
                G_index[*iter2]=1;
            }
            for(iter2=con_Y.begin();iter2!=con_Y.end();++iter2){
                G_index[*iter2]=1;
            }
        }
        if(iteration){
            isMax=false;
        }

    }

    for(iter2=con_Y.begin();iter2!=con_Y.end();++iter2){
        G_index[*iter2]=0;
    }
    for(iter2=cur_Y.begin();iter2!=cur_Y.end();++iter2){
        G_prune[*iter2]=1;
        G_index[*iter2]=0;
    }

    for(iter1=ext_p.begin();iter1!=ext_p.end();++iter1){
        //mark
        ////cout<<*iter1<<endl;
        int temp_count=0;
        bool iteration=true;
        for(iter2=con_Y.begin();iter2!=con_Y.end();++iter2){
            G_index[*iter2]=1;
        }
        for(iter2=cur_Y.begin();iter2!=cur_Y.end();++iter2){
            G_index[*iter2]=1;
        }
        //check new point
        for(int i=0;i<degree[*iter1];++i){
            if(G_index[Graph[*iter1][i]]==1){
                G_index[Graph[*iter1][i]]=2;
                temp_count++;
            }
        }
        if(temp_count>=S-k){
            G_prune[*iter1]=1;
        }
        if(isMax&&temp_count>=temp_num-1-k){
            for(iter2=con_Y.begin();iter2!=con_Y.end();++iter2){
                if(G_index[*iter2]==1&&G_State[*iter2]<(int)X.size()+1-k){
                    iteration=false;
                    break;
                }
            }
            for(iter2=cur_Y.begin();iter2!=cur_Y.end()&&iteration;++iter2){
                if(G_index[*iter2]==1&&G_State[*iter2]<(int)X.size()+1-k){
                    iteration=false;
                    break;
                }
            }
        }else{
            iteration=false;
        }
        if(iteration){
            isMax=false;
            //return false;
        }
    }
    for(iter2=con_Y.begin();iter2!=con_Y.end();++iter2){
        G_index[*iter2]=0;
    }
    for(iter2=cur_Y.begin();iter2!=cur_Y.end();++iter2){
        G_index[*iter2]=0;
    }
    return isMax;

}

void biplexminer::Special(vector<int> &X, vector<int> &var_Y, vector<int> &con_Y, vector<int> &cur_Y,vector<int> &candY_q, vector<int> &ext_p, vector<int> &ext_q){

    vector<int>::iterator iter1;
    vector<int>::iterator iter;
    if(k==1){
        for(iter1=var_Y.begin();iter1!=var_Y.end();++iter1){
            cur_Y.push_back(*iter1);
            ////cout<<"())()"<<" "<<ext_q.size()<<" "<<ext_q.front()<<endl;
            bool res=CheckMaximality(X,con_Y,cur_Y,candY_q,ext_p,ext_q);
            ////cout<<"((("<<endl;
            if(res){
                res_count++;
		    if(res_size<GetSize(X,cur_Y,con_Y))
            {
                res_size = GetSize(X,cur_Y,con_Y);
                Lnode=X;
                Rnode=cur_Y;
                for(int ri=0;ri<con_Y.size();ri++)
                {
                    int curele=con_Y[ri];
                    Rnode.push_back(curele);
                }
            }
	  	        
                //Output(X,cur_Y,con_Y);

            }
            cur_Y.pop_back();
        }
    }else if(k==2){
        vector<int>::iterator iter2;
        for(iter1=var_Y.begin();iter1!=var_Y.end();++iter1){
            iter2=iter1;iter2++;
            cur_Y.push_back(*iter1);
            for(;iter2!=var_Y.end();++iter2){
                cur_Y.push_back(*iter2);
                bool res=CheckMaximality(X,con_Y,cur_Y,candY_q,ext_p,ext_q);
                if(res){
                    res_count++;
		    if(res_size<GetSize(X,cur_Y,con_Y))
            {
                res_size = GetSize(X,cur_Y,con_Y);
                Lnode=X;
                Rnode=cur_Y;
                for(int ri=0;ri<con_Y.size();ri++)
                {
                    int curele=con_Y[ri];
                    Rnode.push_back(curele);
                }
            }
	    		

                    //Output(X,cur_Y,con_Y);
                }
                cur_Y.pop_back();
            }
            cur_Y.pop_back();
        }

    }else if(k==3){
        vector<int>::iterator iter2;
        vector<int>::iterator iter3;
        for(iter1=var_Y.begin();iter1!=var_Y.end();++iter1){
            iter2=iter1;iter2++;
            cur_Y.push_back(*iter1);
            for(;iter2!=var_Y.end();iter2++){
                iter3=iter2;iter3++;
                cur_Y.push_back(*iter2);
                for(;iter3!=var_Y.end();++iter3){
                    cur_Y.push_back(*iter3);
                    bool res=CheckMaximality(X,con_Y,cur_Y,candY_q,ext_p,ext_q);
                    if(res){
                        res_count++;
			if(res_size<GetSize(X,cur_Y,con_Y))
			    {res_size = GetSize(X,cur_Y,con_Y);
                Lnode=X;
                Rnode=cur_Y;
                for(int ri=0;ri<con_Y.size();ri++)
                {
                    int curele=con_Y[ri];
                    Rnode.push_back(curele);
                }}
                        //Output(X,cur_Y,con_Y);

                    }
                    cur_Y.pop_back();
                }
                cur_Y.pop_back();
            }
            cur_Y.pop_back();
        }


    }else if(k==4){
        vector<int>::iterator iter2;
        vector<int>::iterator iter3;
        vector<int>::iterator iter4;

        for(iter1=var_Y.begin();iter1!=var_Y.end();++iter1){
            iter2=iter1;iter2++;
            cur_Y.push_back(*iter1);
            for(;iter2!=var_Y.end();iter2++){
                iter3=iter2;iter3++;
                cur_Y.push_back(*iter2);
                for(;iter3!=var_Y.end();++iter3){
                    iter4=iter3;iter4++;
                    cur_Y.push_back(*iter3);
                    for(;iter4!=var_Y.end();++iter4){
                        cur_Y.push_back(*iter4);
                        bool res=CheckMaximality(X,con_Y,cur_Y,candY_q,ext_p,ext_q);
                        if(res){
                            res_count++;
			    if(res_size<GetSize(X,cur_Y,con_Y)){
			         res_size = GetSize(X,cur_Y,con_Y);
                     Lnode=X;
                Rnode=cur_Y;
                for(int ri=0;ri<con_Y.size();ri++)
                {
                    int curele=con_Y[ri];
                    Rnode.push_back(curele);
                }}
                            //Output(X,cur_Y,con_Y);
                        }
                        cur_Y.pop_back();
                    }
                    cur_Y.pop_back();
                }
                cur_Y.pop_back();
            }
            cur_Y.pop_back();
        }

    }
}

void biplexminer::Output(vector<int> &L, vector<int> &R1, vector<int> &R2){
    vector<int>::iterator out_iter1, out_iter2;
    cout<<"L: ";
    for(out_iter1=L.begin(),out_iter2=L.end();out_iter1!=out_iter2;++out_iter1){
        cout<<*out_iter1<<" ";
    }
    cout<<"R: ";
    for(out_iter1=R1.begin(),out_iter2=R1.end();out_iter1!=out_iter2;++out_iter1){
        cout<<*out_iter1<<" ";
    }
    for(out_iter1=R2.begin(),out_iter2=R2.end();out_iter1!=out_iter2;++out_iter1){
        cout<<*out_iter1<<" ";
    }
    cout<<endl;
}

int biplexminer::GetSize(vector<int> &L, vector<int> &R1, vector<int> &R2){
    // //cout<<"now print current Graph: "<<endl;
    // //cout<<Graph_size<<' ';
    // //cout<<Bipartite_index<<' ';
    // //cout<<10000<<endl;
    // for(int i=0;i<Graph_size;i++)
    // {
    //     //cout<<i;
    //     for(int j=0;j<degree[i];j++)
    //     //cout<<' '<<Graph[i][j];
    //     //cout<<endl;
    // }

    int res = 0;
    int *G_count;
    G_count=new int[Graph_size];
    for(int i=Bipartite_index;i<Graph_size;++i){
        G_count[i] = 0;
    }

    for(auto i:L){
        for(int j=0;j<degree[i];++j){
            G_count[Graph[i][j]]+=1;
        }
    }

    for(auto i:R1)
        res += G_count[i];

    for(auto i:R2)
        res+= G_count[i];
    
    if(res > 10){
        ////cout<<"L: "<<L.size()<<", R: "<<R1.size()<<", R2: "<<R2.size()<<endl;
    vector<int>::iterator out_iter1, out_iter2;
    //cout<<"L: ";
    for(out_iter1=L.begin(),out_iter2=L.end();out_iter1!=out_iter2;++out_iter1){
        //cout<<*out_iter1<<" ";
    }
    //cout<<"R: ";
    for(out_iter1=R1.begin(),out_iter2=R1.end();out_iter1!=out_iter2;++out_iter1){
        //cout<<*out_iter1<<" ";
    }
    for(out_iter1=R2.begin(),out_iter2=R2.end();out_iter1!=out_iter2;++out_iter1){
        //cout<<*out_iter1<<" ";
    }
    //cout<<" size: "<<res;
    //cout<<endl;
    }
    

    return res;
}

