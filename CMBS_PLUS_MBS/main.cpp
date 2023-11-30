#include<iostream>
#include<cstring>
#include<time.h>
#include<list>
#include<vector>
#include "args.hxx"
#include "Util.h"
#include "core.h"
#include "biplex.h"
#include <chrono>
using namespace std;
typedef std::chrono::duration<double> tms;
#define FILELEN 1024
int main(int argc, char** argv) {

	char filepath[1024] = ".........";
	int k = 1;
	args::ArgumentParser parser(
        "CoreBB a algorithm to search the maximum k-biplex in bipartite graph\n");
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
	if (k<=0) {
		fprintf(stderr, "k, theta and num should be at least 1\n");
		exit(-1);
	}

	Util util;
	int bi=0;
    int *degree=NULL;
    int **Graph=NULL;
    cout<<"read graph"<<endl;
    int graph_size=util.ReadGraph(filepath,Graph,degree,bi);
    cout<<"graph is loaded"<<endl;
    Core core(Graph,degree,graph_size,bi,k);
    XYPair XYPairs=core.GetPairs(theta_l,theta_r);
    cout<<"core decomposition"<<endl;
    int maximum=theta_l*(theta_r-k);
    int xmax=0;
    int ymax=0;
    xmax=core.GetMaxX(theta_r);
    ymax=core.GetMaxY(theta_l);
    cout<<"xmax is "<<xmax<<" / "<<XYPairs.GetMaxX()<<" / "<<XYPairs.SearchY(theta_r).upper<<endl;
    cout<<"ymax is "<<ymax<<" / "<<XYPairs.GetMaxY()<<" / "<<XYPairs.SearchX(theta_l).upper<<endl;
    printf("xmax: %d  ymax: %d \n",xmax,ymax);
    bool flag1=false;
    cout<<"now find reduction search core space"<<endl;
    const auto begin=std::chrono::steady_clock().now();
    int upper=ymax;
    int low=(ymax/10)>theta_r?(ymax/10):theta_r;
    // XYPairs.printXYPair();
    cout<<"Y-side reduction: "<<endl;        
    while(true)
    {
        int x=theta_l;
        int y=(upper+low)/2;
        printf("x: %d  y: %d \n",x,y);
        int ** pg;
        int * pd;
        int pb=0;
        const auto core_prune_begin=std::chrono::steady_clock().now();
        int pgs=core.XYCorePrune(pg,pd,pb,y-k,x-k);
        const auto core_prune_end=std::chrono::steady_clock().now();
        tms core_prune_time=core_prune_end-core_prune_begin;
        printf("Core prunned! and cost time %lf \n",core_prune_time.count());
        const auto search_begin=std::chrono::steady_clock().now();
        biplexminer_xy bminer_xy(pg,pgs,pb,pd,x,y,k);
        bminer_xy.Miner();
        const auto search_end=std::chrono::steady_clock().now();
        tms search_time=search_end-search_begin;
        printf("Search finished! and cost time %lf \n",search_time.count());
        if(bminer_xy.res_size!=0)
        {
            if(maximum<bminer_xy.res_size)
            {
                maximum=bminer_xy.res_size;
                XYPairs.update(maximum);
                // XYPairs.printXYPair();
            }
            low=y;
        }
        else
        {
            if((x!=theta_l)||(y!=theta_r))
            {
                XYPairs.removeY(y);
            }
            if(upper>y) upper=y;
        }
        delete [] pd;
        for(int d=0;d<pgs;d++)
        {
            delete [] pg[d];
        }
        delete [] pg;
        if(upper==low||upper-low==1||low>upper) break;
    }
    cout<<"X-side reduction: "<<endl;
    int ymin=XYPairs.GetMinY();
    xmax=XYPairs.GetMaxX();
    cout<<"new xmax is : "<<xmax<<endl;
    flag1=false;
    upper=xmax;
    low=(xmax/10)>theta_l?(xmax/10):theta_l;
    while(true&&xmax>theta_l)
    {
        int x=(upper+low)/2;
        int y=ymin;
        printf("x: %d  y: %d \n",x,y);       
        int ** pg;
        int * pd;
        int pb=0;
        const auto core_prune_begin=std::chrono::steady_clock().now();
        int pgs=core.XYCorePrune(pg,pd,pb,y-k,x-k);
        const auto core_prune_end=std::chrono::steady_clock().now();
        tms core_prune_time=core_prune_end-core_prune_begin;
        printf("Core prunned! and cost time %lf \n",core_prune_time.count());
        const auto search_begin=std::chrono::steady_clock().now();
        biplexminer_xy bminer_xy(pg,pgs,pb,pd,x,y,k);
        bminer_xy.Miner();
        const auto search_end=std::chrono::steady_clock().now();
        tms search_time=search_end-search_begin;
        printf("Search finished! and cost time %lf \n",search_time.count());
        if(bminer_xy.res_size!=0)
        {
            if(maximum<bminer_xy.res_size)
            {
                maximum=bminer_xy.res_size;
                XYPairs.update(maximum);
                // XYPairs.printXYPair();
            }
            low=x;
        }
        else
        {
            if((x!=theta_l)||(y!=theta_r))
            {
                XYPairs.removeX(x);
            }
            if(upper>x) upper=x;
            
        }
        delete [] pd;
        for(int d=0;d<pgs;d++)
        {
            delete [] pg[d];
        }
        delete [] pg;
        if(upper==low||upper-low==1||low>upper) break;
    }
    const auto end=std::chrono::steady_clock().now();
    tms total_time=end-begin;
    printf("Total time: %lf \n",total_time.count());
    cout<<"now print the maximum : "<<maximum<<endl;
	return 0;
}