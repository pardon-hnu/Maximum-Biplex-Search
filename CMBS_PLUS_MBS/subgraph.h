struct subgraph{
    int edges;
    int left_size, right_size;
    vector<int> left; vector<int> right;
    subgraph(){
        this->edges=0;
        this->left_size=0;
        this->right_size=0;
    }
    subgraph(int edges, int left_size, int right_size){
        this->edges=edges;
        this->left_size=left_size;
        this->right_size=right_size;
    }
    bool operator>(const subgraph& obj) const{
        return this->edges>obj.edges;
    }
    void Sort(){
        sort(left.begin(),left.end());
        sort(right.begin(),right.end());
    }
    friend ostream &operator<<( ostream &output, const subgraph& obj ){ 
        output<<"|E|: "<<obj.edges<<"  |L|: "<<obj.left_size<<"  |R|: "<<obj.right_size<<endl;
        output<<"L: ";
        for(int i:obj.left) output<<i<<" ";
        output<<endl;
        output<<"R: ";
        for(int i:obj.right) output<<i<<" ";
        output<<endl;
        return output;
    }
};