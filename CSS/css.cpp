#include<iostream>
#include<string>
#include<time.h>
#include<list>
#include<vector>
#include<cstring>
#include "args.hxx"
#include"Util.h"
#include"Framework.h"
#include"Corepruning.h"

#define FILELEN 1024


int main(int argc, char** argv) {

	char filepath[1024] = ".........";

	int k = 1;

	args::ArgumentParser parser(
        "iTraversal, an algorithm for enumerating all maximal biplexes\n");

    args::HelpFlag help(parser, "help", "Display this help menu",
                        {'h', "help"});
    args::Group required(parser, "", args::Group::Validators::All);

    args::ValueFlag<std::string> benchmark_file(
        parser, "benchmark", "Path to benchmark", {'f', "file"},
        "");

    args::ValueFlag<int> Threshold_l(parser, "Lower Bound", "The lower bound of the size of k-biplex", {'u', "u"}, 1);

    args::ValueFlag<int> Threshold_r(parser, "Lower Bound", "The lower bound of the size of k-biplex", {'v', "v"}, 1);

	args::ValueFlag<int> K(parser, "para k", "The parameter k", {'k', "k"}, 1);


    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 0;
    } catch (args::ValidationError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 0;
    }

	strncpy(filepath, args::get(benchmark_file).c_str(), FILELEN);
	k = args::get(K);
    int theta_l = args::get(Threshold_l);
    int theta_r = args::get(Threshold_r);
    if(theta_l<2*k+1){
        cout<<"left size should be larger than "<<2*k+1<<endl;
        exit(-1);
    }
    if(theta_r<2*k+1){
        cout<<"right size should be larger than "<<2*k+1<<endl;
        exit(-1);
    }
    int top_k=1;
    if(top_k<=0){
        cout<<"K should be a positive integer"<<endl;
        exit(-1); 
    }

	if (k<0) {
		fprintf(stderr, "k, theta and num should be at least 1\n");
		exit(-1);
	}


	Util util;
	int bi=0;
    int *degree=NULL;
    int **Graph=NULL;
    int graph_size=util.ReadGraph(filepath,Graph,degree,bi);
	
    CoreLocate core(Graph, degree, graph_size,bi,2*k);
    core.Coredecompose();
    int **pG=NULL;
    int *pd=NULL;
    int pb=0;
    int pgs=core.CorePrune(pG,pd,pb);
    
    Framework bkmb(pG,pd,pgs,pb,k,theta_l,theta_r,top_k);
    time_t s1 = clock();
    bkmb.test(theta_l,theta_r);
    time_t s2 = clock();
    cout<<"Running Time: "<<((double)(s2-s1)/CLOCKS_PER_SEC)<<" sec"<<endl;
    // cout<<bkmb.MaxBP<<endl;
    while(!bkmb.res_top_k.empty()){
        subgraph g=bkmb.res_top_k.top();
        cout<<g<<endl;
        bkmb.res_top_k.pop();
    }
	return 0;
}