#include "ap_int.h"
#include "hls_vector.h"
#include "hls_stream.h"
#include "hls_math.h"
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <chrono>
using namespace std;
using namespace hls;

//control define

#define CONTROL
//#define PRINT
//#define DEBUG_CORE
//#define DEBUG_REDUCTION
//const define
//#define INFO


#define MAX 			64
#define GRAPH_SIZE 		44
#define EDGE_SIZE		99
#define MAX_DEGREE		12
#define MAX_TRIPCOUNT	1000000
#define BATCH_SIZE		10
#define DDR_SIZE	 	100
#define BRAM_SIZE		100

#define NODE_WIDTH		32
#define DEGREE_WIDTH	16
//type define
typedef ap_int<NODE_WIDTH> nodeT;
typedef ap_int<DEGREE_WIDTH> degreeT;
typedef hls::vector<nodeT,MAX> set_node;
struct parvector
{
	set_node data;
	nodeT size;
	parvector()
	{

	}
	parvector(set_node d,nodeT s)
	{
		data=d;
		size=s;
	}
	void setvalue(parvector p)
	{
		data=p.data;
		size=p.size;
	}
	void init()
	{
#pragma HLS PIPELINE
		size=0;
	}
	void add(nodeT v,nodeT * pos)
	{
		data[size]=v;
		pos[v]=size;
		size++;
	}
	void insert(nodeT v)
	{
		data[size]=v;
		size++;
	}
	void remove(nodeT v,nodeT * pos)
	{
		nodeT v_pos=pos[v];
		nodeT last_id=data[size-1];
		data[v_pos]=last_id;
		pos[last_id]=v_pos;
		pos[v]=GRAPH_SIZE;
		size--;
	}
	bool exist(nodeT v,nodeT * pos)
	{
		bool flag;
		if(pos[v]==GRAPH_SIZE) flag=false;
		else flag=true;
		return flag;
	}
};

struct branch
{
	parvector P_L;
	parvector P_R;
	parvector C_L;
	parvector C_R;
	parvector E_L;
	parvector E_R;
	branch()
	{

	}
	branch(parvector pl,parvector pr,parvector cl,parvector cr,parvector el,parvector er)
	{
		P_L=pl;
		P_R=pr;
		C_L=cl;
		C_R=cr;
		E_L=el;
		E_R=er;
	}
	void setvalue(parvector pl,parvector pr,parvector cl,parvector cr,parvector el,parvector er)
	{
		P_L=pl;
		P_R=pr;
		C_L=cl;
		C_R=cr;
		E_L=el;
		E_R=er;
	}
	void init()
	{
		P_L.init();
		P_R.init();
		C_L.init();
		C_R.init();
		E_L.init();
		E_R.init();
	}
};
extern "C"
{
void Framework(branch * ddr,nodeT * Index,nodeT * Neighbor,degreeT * Degree,nodeT Graph_size,nodeT Bipartite,nodeT K,nodeT Alpha,nodeT Beta,branch * mbp,nodeT * maximum);
}
