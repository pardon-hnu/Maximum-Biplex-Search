#include "biplex1.h"

biplexminer_xy::biplexminer_xy(int **Graph, int Graph_size, int Bipartite_index,int *degree, int alpha, int beta, int k){
    this->Graph=Graph;
    this->Bipartite_index=Bipartite_index;
    this->alpha=alpha;
    this->beta=beta;
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

biplexminer_xy::~biplexminer_xy(){
    delete []G_State;
    delete []G_index;
    delete []G_label;
    delete []G_temp;
    delete []G_prune;
}

void biplexminer_xy::Miner(){
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
    Enumerate(X,N_X,Cand_exts,cand_q);
}

void biplexminer_xy::Miner1(){
    vector<int> X;
    vector<int> N_X;
    vector<int> cand_q;
    vector<int> Cand_exts;
    for(int i=Bipartite_index;i<Graph_size;++i){
        Cand_exts.push_back(i);
    }
    for (int i = 0; i < Bipartite_index; i++){
        N_X.push_back(i);
    }

    for(int i=Bipartite_index;i<Graph_size;++i){
        G_label[i]=0;
    }

    for(int i=0;i<Bipartite_index;++i){
        G_label[i]=1;
    }


    Enumerate1(X,N_X,Cand_exts,cand_q);
}

void biplexminer_xy::Enumerate(vector<int> &X, vector<int> &N_X, vector<int> &Cand_p, vector<int> &Cand_q){
    int temp_node;
    int temp_size;
    vector<int> new_N_X;
    vector<int> new_X;
    vector<int> new_Cand_p;
    vector<int> new_Cand_q;
    vector<int>:: iterator iter1;
    vector<int>:: iterator iter2;
    //N_X: V侧顶点，Cand_p:U侧顶点

    for(iter1=Cand_q.begin();iter1!=Cand_q.end();++iter1){
        new_Cand_q.push_back(*iter1);
    }
    //copy cand_q

    while(!Cand_p.empty()){

        new_N_X.clear();
        new_X.clear();
        new_Cand_p.clear();


        for(iter1=X.begin();iter1!=X.end();++iter1){
            new_X.push_back(*iter1);
        }

        //ȡ��cand_p�еĽڵ㲢�������new_x��
        int temp_node=Cand_p.front();
        new_X.push_back(temp_node);
        //?
        G_label[temp_node]=1;
        Cand_p.erase(Cand_p.begin());

        if((int)(new_X.size())<= k+1&&(int)(new_X.size())<alpha){
            //����¼���ڵ���ھӵ�state, ����ֵ��ʾ��X�нڵ������ĸ���
            //update G_State
            for(int i=0;i<degree[temp_node];++i){
                G_State[Graph[temp_node][i]]+=1;
            }

            //ͨ��state����X���ھӼ���new_N_X
            //Create new_N_X //improve
            if(new_X.size()== k+1){
                for(int i=Bipartite_index;i<Graph_size;++i){
                    if(G_State[i]>0)
                        new_N_X.push_back(i);
                    else
                        G_label[i]=0;
                }
            }else{
                new_N_X=N_X;
            }

            //update G_State
            temp_size=degree[temp_node];
            G_State[temp_node]=temp_size;


            //create new cand_p
            for(iter1=Cand_p.begin();iter1!=Cand_p.end();++iter1){

                //cand�еĽڵ������������
                temp_size=degree[*iter1];
                if(temp_size<beta-k){
                    continue;
                }

                //possible debug label
                for(int i=0;i<temp_size;++i){
                    if(G_label[Graph[*iter1][i]]==1)
                        G_index[Graph[*iter1][i]]=1;
                }
                int temp_count=0;
                for(int i=0;i<degree[temp_node];++i){
                    if(G_index[Graph[temp_node][i]]==1){
                        temp_count++;
                    }
                }

                if(temp_count<beta-2*k){
                    for(int i=0;i<temp_size;++i){
                        G_index[Graph[*iter1][i]]=0;
                    }
                    continue;
                }
                for(int i=0;i<temp_size;++i){
                    G_index[Graph[*iter1][i]]=0;
                }
                new_Cand_p.push_back(*iter1);

            }

            //check the probability of extension
            if(degree[temp_node]>=beta-k){
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
                    if(temp_count<beta-2*k){
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
        }else{
            //update degree of Y
            bool cost=false;
            int temp_count=0;
            for(int i=0;i<degree[temp_node];++i){
                if(G_label[Graph[temp_node][i]]==1){
                    temp_count++;
                    G_State[Graph[temp_node][i]]+=1;
                }
            }
            //update degree of X
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
            //pruning 1 : check probability of extension
            if(new_N_X.size()<0){
                Extension=false;
            }

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
                        if(temp_count<beta-2*k){
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
            ////cout<<"prunning2 : "<<Extension<<endl;
            //pruning 3
            vector<int> Component;
            if(Extension){

                for(iter1=Cand_p.begin();iter1!=Cand_p.end();++iter1){
                    if(degree[*iter1]>=beta-k){
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
                        }else if(temp_count>=beta-k){
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
                        if(temp_count<beta-2*k){
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

                if(new_X.size()>=alpha){
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

                if(new_X.size()+new_Cand_p.size()>=alpha&&new_N_X.size()>=alpha){
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

void biplexminer_xy::Enumerate1(vector<int> &X, vector<int> &N_X, vector<int> &Cand_p, vector<int> &Cand_q){
    int temp_node;
    int temp_size;
    vector<int> new_N_X;
    vector<int> new_X;
    vector<int> new_Cand_p;
    vector<int> new_Cand_q;
    vector<int>:: iterator iter1;
    vector<int>:: iterator iter2;



    for(iter1=Cand_q.begin();iter1!=Cand_q.end();++iter1){
        new_Cand_q.push_back(*iter1);
    }

    while(!Cand_p.empty()){

        new_N_X.clear();
        new_X.clear();
        new_Cand_p.clear();


        for(iter1=X.begin();iter1!=X.end();++iter1){
            new_X.push_back(*iter1);
        }

        //ȡ��cand_p�еĽڵ㲢�������new_x��
        int temp_node=Cand_p.front();
        new_X.push_back(temp_node);
        //?
        G_label[temp_node]=1;
        Cand_p.erase(Cand_p.begin());

        if((int)(new_X.size())<= k+1&&(int)(new_X.size())<alpha){
            //����¼���ڵ���ھӵ�state, ����ֵ��ʾ��X�нڵ������ĸ���
            //update G_State
            for(int i=0;i<degree[temp_node];++i){
                G_State[Graph[temp_node][i]]+=1;
            }

            //ͨ��state����X���ھӼ���new_N_X
            //Create new_N_X //improve
            if(new_X.size()== k+1){
                for(int i=0;i<Bipartite_index;++i){
                    if(G_State[i]>0)
                        new_N_X.push_back(i);
                    else
                        G_label[i]=0;
                }
            }else{
                new_N_X=N_X;
            }

            //update G_State
            temp_size=degree[temp_node];
            G_State[temp_node]=temp_size;


            //create new cand_p
            for(iter1=Cand_p.begin();iter1!=Cand_p.end();++iter1){

                //cand�еĽڵ������������
                temp_size=degree[*iter1];
                if(temp_size<beta-k){
                    continue;
                }

                //possible debug label
                for(int i=0;i<temp_size;++i){
                    if(G_label[Graph[*iter1][i]]==1)
                        G_index[Graph[*iter1][i]]=1;
                }
                int temp_count=0;
                for(int i=0;i<degree[temp_node];++i){
                    if(G_index[Graph[temp_node][i]]==1){
                        temp_count++;
                    }
                }

                if(temp_count<beta-2*k){
                    for(int i=0;i<temp_size;++i){
                        G_index[Graph[*iter1][i]]=0;
                    }
                    continue;
                }
                for(int i=0;i<temp_size;++i){
                    G_index[Graph[*iter1][i]]=0;
                }
                new_Cand_p.push_back(*iter1);

            }

            //check the probability of extension
            if(degree[temp_node]>=beta-k){
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
                    if(temp_count<beta-2*k){
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
                for(int i=0;i<Bipartite_index;++i){
                    G_label[i]=1;
                }
            }
        }else{
            //update degree of Y
            bool cost=false;
            int temp_count=0;
            for(int i=0;i<degree[temp_node];++i){
                if(G_label[Graph[temp_node][i]]==1){
                    temp_count++;
                    G_State[Graph[temp_node][i]]+=1;
                }
            }
            //update degree of X
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
            //pruning 1 : check probability of extension
            if(new_N_X.size()<0){
                Extension=false;
            }

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
                        if(temp_count<beta-2*k){
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
            ////cout<<"prunning2 : "<<Extension<<endl;
            //pruning 3
            vector<int> Component;
            if(Extension){

                for(iter1=Cand_p.begin();iter1!=Cand_p.end();++iter1){
                    if(degree[*iter1]>=beta-k){
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
                        }else if(temp_count>=beta-k){
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
                        if(temp_count<beta-2*k){
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

                if(new_X.size()>=alpha){
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

                if(new_X.size()+new_Cand_p.size()>=alpha&&new_N_X.size()>=alpha){
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

void biplexminer_xy::Generate(vector<int> &X, vector<int> &N_X, vector<int> &Cand_p, vector<int> &Cand_q){
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
        if(temp_count>=alpha-k){
            ext_q.push_back(*iter1);
        }
    }
    ////cout<<"%%%"<<constant_Y.size()<<" "<<var_Y.size()<<endl;
    if(invalid_X.size()==0&&(constant_Y.size()+var_Y.size()>=beta)){
        bool res=CheckMaximality(X,constant_Y,var_Y,candY_p,Cand_p,ext_q);
        if(res){
            res_count++;
            int temp_size = GetSize(X,var_Y,constant_Y);
            if(res_size < temp_size)
               { res_size = temp_size; Output(X,var_Y,constant_Y);}

            //if(OutputResults)
                //OutputFormat(X,var_Y,constant_Y);
        }
    }else if(invalid_X.size()==1&&k<=4){
         if((int)constant_Y.size()>=beta-k&&(constant_Y.size()+var_Y.size())>=beta){
             Special(X,var_Y,constant_Y,cur_Y,candY_p,Cand_p,ext_q);
         }
    }else{
        List(X,var_Y,constant_Y,cur_Y,candY_p,Cand_p,ext_q);
    }
    //recovery
    for(iter1=X.begin();iter1!=X.end();++iter1){
        G_temp[*iter1]=0;
    }
}


void biplexminer_xy::List(vector<int> &X, vector<int> var_Y, vector<int> &con_Y, vector<int> &cur_Y, vector<int> &candY_q, vector<int> &ext_p,
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
    if(var_Y.empty()&&check&&(cur_Y.size()+con_Y.size())>=beta){
        bool res=CheckMaximality(X,con_Y,cur_Y,candY_q,ext_p,ext_q);
        if(res){
            res_count++;
            int temp_size = GetSize(X,cur_Y,con_Y);
            if(res_size < temp_size)
                {res_size = temp_size;Output(X,cur_Y,con_Y);}

            //if(OutputResults)
                //OutputFormat(X,cur_Y,con_Y);

        }
    }

}

bool biplexminer_xy::CheckMaximality(vector<int> &X, vector<int> &con_Y, vector<int> &cur_Y,
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
        if(temp_count>=alpha-k){
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

void biplexminer_xy::Special(vector<int> &X, vector<int> &var_Y, vector<int> &con_Y, vector<int> &cur_Y,vector<int> &candY_q, vector<int> &ext_p, vector<int> &ext_q){

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
                int temp_size = GetSize(X,cur_Y,con_Y);
                if(res_size < temp_size)
                    {res_size = temp_size; Output(X,cur_Y,con_Y);}

                //if(OutputResults)
                    //OutputFormat(X,cur_Y,con_Y);

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
                    int temp_size = GetSize(X,cur_Y,con_Y);
                    if(res_size < temp_size)
                       { res_size = temp_size; Output(X,cur_Y,con_Y);}

                    //if(OutputResults)
                        //OutputFormat(X,cur_Y,con_Y);
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
                        int temp_size = GetSize(X,cur_Y,con_Y);
                        if(res_size < temp_size)
                            {res_size = temp_size;Output(X,cur_Y,con_Y);}

                        //if(OutputResults)
                            //OutputFormat(X,cur_Y,con_Y);

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
                            int temp_size = GetSize(X,cur_Y,con_Y);
                            if(res_size < temp_size)
                               { res_size = temp_size;Output(X,cur_Y,con_Y);}

                            //if(OutputResults)
                                //OutputFormat(X,cur_Y,con_Y);
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

int biplexminer_xy::GetSize(vector<int> &L, vector<int> &R1, vector<int> &R2){
    int res = 0;
    int *G_count;
    G_count=new int[Graph_size];
    for(int i=Bipartite_index;i<Graph_size;++i){
        G_count[i] = 0;
    }

    vector<int>::iterator out_iter1, out_iter2;
    ////cout<<"L: "<<L.size()<<", R: "<<R1.size()+R2.size()<<endl;
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

    for(auto i:L){
        for(int j=0;j<degree[i];++j){
            G_count[Graph[i][j]]+=1;
        }
    }

    for(auto i:R1)
        res += G_count[i];

    for(auto i:R2)
        res+= G_count[i];

    return res;

    
}

void biplexminer_xy::Output(vector<int> &L, vector<int> &R1, vector<int> &R2){
    vector<int>::iterator out_iter1, out_iter2;
    vector<int> a,b;
    cout<<"L: ";
    for(out_iter1=L.begin(),out_iter2=L.end();out_iter1!=out_iter2;++out_iter1){
        cout<<*out_iter1<<" ";
        a.push_back(*out_iter1);
    }
    cout<<"R: ";
    for(out_iter1=R1.begin(),out_iter2=R1.end();out_iter1!=out_iter2;++out_iter1){
        cout<<*out_iter1<<" ";
        b.push_back(*out_iter1);
    }
    for(out_iter1=R2.begin(),out_iter2=R2.end();out_iter1!=out_iter2;++out_iter1){
        cout<<*out_iter1<<" ";
        b.push_back(*out_iter1);
    }
    maxL=a;
    maxR=b;
    cout<<endl;
}