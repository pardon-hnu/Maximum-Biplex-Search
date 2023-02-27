#include "utils.h"

void coredecompositon(nodeT * index,nodeT * neighbor,degreeT * degree,nodeT graph_size,nodeT bipartite,nodeT * pg,nodeT * pi,nodeT * pn,degreeT * pd,nodeT & pgs,nodeT & pb,nodeT alpha,nodeT beta)
{
#pragma HLS INLINE off
#ifdef __SYNTHESIS__
	nodeT temp_index[GRAPH_SIZE];
	degreeT temp_deg[GRAPH_SIZE];
#else
	nodeT * temp_index=new nodeT[GRAPH_SIZE];
	degreeT * temp_deg=new degreeT[GRAPH_SIZE];
#endif
	Core_INIT:
	for(nodeT i=0;i<GRAPH_SIZE;i++)
	{
#pragma HLS PIPELINE
		temp_index[i]=1;
		temp_deg[i]=degree[i];
	}

	bool flag=true;
	nodeT countL=0, countR=0;
	CORE_DECOMPOSITON:
	while(flag)
	{
		flag=false;
		CORE_ITERATION:
		for(nodeT temp_i=0;temp_i<GRAPH_SIZE;temp_i++)
		{
#pragma HLS PIPELINE
			if(temp_i<bipartite)
			{
				if(temp_index[temp_i]&&temp_deg[temp_i]<alpha)
				{
					temp_index[temp_i]=0;
					flag=true;
					countL++;
					nodeT base=index[temp_i];
					CORE_DELETE_L:
					for(nodeT temp_j=0;temp_j<MAX_DEGREE;temp_j++)
					{
						if(temp_j<degree[temp_i])
						{
							nodeT offset=base+temp_j;
							nodeT temp_node3=neighbor[offset];
							temp_deg[temp_node3]--;
						}
					}
				}
			}
			else if(temp_i>=bipartite&&temp_i<graph_size)
			{
				if(temp_index[temp_i]&&temp_deg[temp_i]<beta)
				{
					temp_index[temp_i]=0;
					flag=true;
					countR++;
					nodeT base=index[temp_i];
					CORE_DELETE_R:
					for(nodeT temp_j=0;temp_j<MAX_DEGREE;temp_j++)
					{
						if(temp_j<degree[temp_i])
						{
							nodeT offset=base+temp_j;
							nodeT temp_node3=neighbor[offset];
							temp_deg[temp_node3]--;
						}
					}
				}
			}
		}
	}
	pgs=graph_size-countL-countR;
	pb=bipartite-countL;
	nodeT temp=0;
	nodeT temp_node=0;
	GIVE_INDEX:
	for(nodeT temp_i=0;temp_i<GRAPH_SIZE;temp_i++)
	{
#pragma HLS PIPELINE
		if(temp_index[temp_i]&&temp_i<graph_size)
		{
			pi[temp_node]=temp;
			pg[temp_node]=temp_i;
			temp_node++;
			temp_index[temp_i]=temp_node;
			temp+=temp_deg[temp_i];
		}
	}
	pi[temp_node]=temp;
	temp=0;
	GIVE_DEGREE:
	for(nodeT temp_i=0;temp_i<GRAPH_SIZE;temp_i++)
	{
#pragma HLS PIPELINE
		if(temp_i<graph_size)
		{
			temp_node=temp_index[temp_i]-1;
			if(temp_node>=0)
			{
				pd[temp_node]=temp_deg[temp_i];
				nodeT base=index[temp_i];
				GIVE_NEIGHBOR:
				for(nodeT temp_j=0;temp_j<MAX_DEGREE;temp_j++)
				{
					if(temp_j<degree[temp_i])
					{
						nodeT temp_node2=neighbor[base+temp_j];

						if(temp_index[temp_node2])
						{
							pn[temp]=temp_index[temp_node2]-1;
							temp++;
						}
					}
				}
			}
		}
	}

}
void init_deal_space(degreeT * degInP,degreeT * degInG,nodeT * G_index,nodeT * Pos1,nodeT * Pos2,nodeT * Pos3,nodeT * Pos4,nodeT * Pos5,nodeT * Pos6)
{
	INIT_DEGINP:
	for(int init_i=0;init_i<GRAPH_SIZE;init_i++)
	{
#pragma HLS PIPELINE
		degInP[init_i]=0;
	}
	INIT_DEGING:
	for(int init_i=0;init_i<GRAPH_SIZE;init_i++)
	{
#pragma HLS PIPELINE
		degInG[init_i]=0;
	}
	INIT_INDEX_G:
	for(int init_i=0;init_i<GRAPH_SIZE;init_i++)
	{
#pragma HLS PIPELINE
		G_index[init_i]=0;
	}
	INIT_POS_1:
	for(int init_i=0;init_i<GRAPH_SIZE;init_i++)
	{
#pragma HLS PIPELINE
		Pos1[init_i]=GRAPH_SIZE;
	}
	INIT_POS_2:
	for(int init_i=0;init_i<GRAPH_SIZE;init_i++)
	{
#pragma HLS PIPELINE
		Pos2[init_i]=GRAPH_SIZE;
	}
	INIT_POS_3:
	for(int init_i=0;init_i<GRAPH_SIZE;init_i++)
	{
#pragma HLS PIPELINE
		Pos3[init_i]=GRAPH_SIZE;
	}
	INIT_POS_4:
	for(int init_i=0;init_i<GRAPH_SIZE;init_i++)
	{
#pragma HLS PIPELINE
		Pos4[init_i]=GRAPH_SIZE;
	}
	INIT_POS_5:
	for(int init_i=0;init_i<GRAPH_SIZE;init_i++)
	{
#pragma HLS PIPELINE
		Pos5[init_i]=GRAPH_SIZE;
	}
	INIT_POS_6:
	for(int init_i=0;init_i<GRAPH_SIZE;init_i++)
	{
#pragma HLS PIPELINE
		Pos6[init_i]=GRAPH_SIZE;
	}
}
void loadpos(parvector s,nodeT * pos)
{
	nodeT temp;
	for(nodeT i=0;i<MAX;i++)
	{
		if(i<s.size)
		{
			temp=s.data[i];
			pos[temp]=i;
		}
	}
}
void Write_bram_to_ddr(branch * ddr,branch * bram,nodeT offset,nodeT lenth)
{
	memcpy((void*)&ddr[offset],(void*)bram,lenth*sizeof(branch));
}
void Read_ddr_to_bram(branch * ddr,branch * bram,nodeT offset,nodeT lenth)
{
	memcpy((void*)bram,(void*)&ddr[offset],lenth*sizeof(branch));
}
void store_branch(branch b,branch * ddr,nodeT & ddr_size,branch * bram,nodeT & bram_size)
{
#pragma HLS INLINE off
	if(bram_size>=BRAM_SIZE)
	{
		Write_bram_to_ddr(ddr,bram,ddr_size,BRAM_SIZE);
		bram_size-=BRAM_SIZE;
		ddr_size+=BRAM_SIZE;
	}
	bram[bram_size]=b;
	bram_size++;
}
void get_branch(branch * batch,nodeT & batch_size,branch * bram,nodeT & bram_size,branch * ddr,nodeT & ddr_size)
{
#pragma HLS INLINE off
	if(bram_size==0)
	{
		nodeT offset;
		nodeT lenth;
		if(ddr_size<BRAM_SIZE)
		{
			offset=0;
			lenth=ddr_size;
		}
		else
		{
			offset=ddr_size-BRAM_SIZE;
			lenth=BRAM_SIZE;
		}
		Read_ddr_to_bram(ddr,bram,offset,lenth);
		ddr_size-=lenth;
		bram_size+=lenth;
	}
	batch_size=0;
	LOAD_BRAM_BRANCH:
	for(nodeT i=0;i<BATCH_SIZE;i++)
	{
#pragma HLS PIPELINE
		if(i<bram_size)
		{
			batch[i]=bram[bram_size-i-1];
			batch_size++;
		}
	}
	bram_size-=batch_size;
}
bool verification(branch b,degreeT * degInG,nodeT left,nodeT right,nodeT last,branch & mbp,nodeT & maximum)
{
#pragma HLS INLINE off
	if(b.P_L.size>=last) return false;
	if(b.P_L.size+b.C_L.size<left||b.P_R.size+b.C_R.size<right) return false;
	if(((b.P_L.size+b.C_L.size)*(b.P_R.size+b.C_R.size)<=maximum)&&((b.P_L.size+b.C_L.size)*(b.P_R.size+b.C_R.size)>0)) return false;
	degreeT Pdegree=0;
	degreeT Gdegree=0;
	nodeT temp_node;
	nodeT temp_i;
	SUM_DEGREE_P:
	for(temp_i=MAX-1;temp_i>=0;temp_i--)
	{
		if(temp_i<b.P_L.size)
		{
			temp_node=b.P_L.data[temp_i];
			Pdegree+=degInG[temp_node];
		}
	}
	Gdegree=Pdegree;
	SUM_DEGREE_C:
	for(temp_i=MAX-1;temp_i>=0;temp_i--)
	{
		if(temp_i<b.C_L.size)
		{
			temp_node=b.C_L.data[temp_i];
			Gdegree+=degInG[temp_node];
		}
	}
	if(Gdegree<=maximum) return false;
	if(b.C_R.size<=0&&b.C_L.size<=0)
	{
		if(b.E_L.size<=0&&b.E_R.size<=0)
		{
			if(Pdegree>maximum)
			{
				mbp=b;
				maximum=Pdegree;
			}
		}
		return false;
	}
	MAXIMAL_CHECK_L:
	for(temp_i=MAX-1;temp_i>=0;temp_i--)
	{
		if(temp_i<b.E_L.size)
		{
			temp_node=b.E_L.data[temp_i];
			if(degInG[temp_node]==b.P_R.size+b.C_R.size)
			{
				return false;
			}
		}
	}
	MAXIMAL_CHECK_R:
	for(temp_i=MAX-1;temp_i>=0;temp_i--)
	{
		if(temp_i<b.E_R.size)
		{
			temp_node=b.E_R.data[temp_i];
			if(degInG[temp_node]==b.P_L.size+b.C_L.size)
			{
				return false;
			}
		}
	}
	return true;
}
bool bit_G[GRAPH_SIZE][GRAPH_SIZE];
void expand1(branch * ddr,nodeT & ddrsize,branch * buffer_in_bram,nodeT & bramsize,
		degreeT * degInG,degreeT * degInP,nodeT * G_index,
		nodeT * pi,nodeT *pn,degreeT * pd,nodeT pgs,nodeT pb,nodeT left,nodeT right,nodeT k,nodeT last,branch & mbp,nodeT & maximum,
		nodeT pivot,parvector P1,parvector P2,parvector C1,parvector C2,parvector E1,parvector E2,nodeT temp_stand1,nodeT temp_stand2,nodeT * PosP1,nodeT * PosP2,nodeT * PosC1,nodeT * PosC2,nodeT * PosE1,nodeT * PosE2)
{
#pragma HLS INLINE off
	nodeT temp_i;
	nodeT temp_j;
	nodeT temp_k;
	nodeT count;
	nodeT temp_value;
	nodeT temp_size1;
	nodeT temp_size2;
	nodeT temp_node;
	nodeT temp_node2;
	nodeT temp_node3;
	nodeT temp_nodei;
	parvector temp_vector;
	temp_vector.init();
	nodeT offset;
	branch new_temp;
	parvector doing;
	doing.init();
	PIVOT_LOOP_1:
	for(temp_i=MAX-1;temp_i>=0;temp_i--)
	{
#pragma HLS PIPELINE
		if(temp_i<C2.size)
		{
			temp_node=C2.data[temp_i];
			if(!bit_G[pivot][temp_node])
			{
				doing.insert(temp_node);
			}
		}
	}
	nodeT p=k-(P2.size-degInP[pivot]);
	nodeT idx=1;
	nodeT rec_a,rec_b;
	rec_a=doing.data[0];
	C2.remove(rec_a,PosC2);
	E2.add(rec_a,PosE2);
	PIVOT_LOOP_2:
	for(temp_i=0;temp_i<MAX_DEGREE;temp_i++)
	{
#pragma HLS PIPELINE
		if(temp_i<pd[rec_a])
		{
			offset=pi[rec_a]+temp_i;
			temp_node2=pn[offset];
			degInG[temp_node2]--;
		}
	}
	if(pivot<pb) new_temp.setvalue(P1,P2,C1,C2,E1,E2);
	else new_temp.setvalue(P2,P1,C2,C1,E2,E1);
	if(verification(new_temp,degInG,left,right,last,mbp,maximum))
	{store_branch(new_temp,ddr,ddrsize,buffer_in_bram,bramsize);}
	E2.remove(rec_a,PosE2);
	PIVOT_LOOP_3:
	while(idx<p+1)
	{
		P2.add(rec_a,PosP2);
		PIVOT_3_1:
		for(temp_i=0;temp_i<MAX_DEGREE;temp_i++)
		{
#pragma HLS PIPELINE
			if(temp_i<pd[rec_a])
			{
				offset=pi[rec_a]+temp_i;
				temp_node=pn[offset];
				degInG[temp_node]++;
				degInP[temp_node]++;
			}
		}
		if(degInP[rec_a]==P1.size-k)
		{
			PIVOT_3_2:
			for(temp_i=MAX-1;temp_i>=0;temp_i--)
			{
#pragma HLS PIPELINE
				if(temp_i<C1.size)
				{
					temp_node=C1.data[temp_i];
					if(!bit_G[rec_a][temp_node]||degInG[temp_node]<temp_stand2-k)
					{
						C1.remove(temp_node,PosC1);
						PIVOT_3_2_1:
						for(temp_j=0;temp_j<MAX_DEGREE;temp_j++)
						{
							if(temp_j<pd[temp_node])
							{
								offset=pi[temp_node]+temp_j;
								temp_node2=pn[offset];
								degInG[temp_node2]--;
							}
						}
					}
				}
			}
			PIVOT_3_3:
			for(temp_i=MAX-1;temp_i>=0;temp_i--)
			{
#pragma HLS PIPELINE
				if(temp_i<E1.size)
				{
					temp_node=E1.data[temp_i];
					if(!bit_G[rec_a][temp_node]||degInG[temp_node]<temp_stand2-k)
					{
						E1.remove(temp_node,PosE1);
					}
				}
			}
		}
		else
		{
			temp_size1=P2.size-k;
			PIVOT_3_4:
			for(temp_i=MAX-1;temp_i>=0;temp_i--)
			{
#pragma HLS PIPELINE
				if(temp_i<C1.size)
				{
					temp_node=C1.data[temp_i];
					if(degInP[temp_node]<temp_size1||degInG[temp_node]<temp_stand2-k)
					{
						C1.remove(temp_node,PosC1);
						PIVOT_3_4_1:
						for(temp_j=0;temp_j<MAX_DEGREE;temp_j++)
						{
							if(temp_j<pd[temp_node])
							{
								offset=pi[temp_node]+temp_j;
								temp_node2=pn[offset];
								degInG[temp_node2]--;
							}
						}
					}
				}
			}
			PIVOT_3_5:
			for(temp_i=MAX-1;temp_i>=0;temp_i--)
			{
#pragma HLS PIPELINE
				if(temp_i<E1.size)
				{
					temp_node=E1.data[temp_i];
					if(degInP[temp_node]<temp_size1||degInG[temp_node]<temp_stand2-k)
					{
						E1.remove(temp_node,PosE1);
					}
				}
			}
		}
		if(C1.size+P1.size<temp_stand1) break;
		temp_vector.init();
		temp_size1=P2.size;
		temp_size2=P1.size;
		PIVOT_3_6:
		for(temp_i=MAX-1;temp_i>=0;temp_i--)
		{
#pragma HLS PIPELINE
			if(temp_i<P1.size)
			{
				temp_node=P1.data[temp_i];
				if(degInP[temp_node]==temp_size1-k&&!bit_G[rec_a][temp_node])
				{
					temp_vector.insert(temp_node);
				}
			}
		}
		PIVOT_3_7:
		for(temp_i=0;temp_i<MAX_DEGREE;temp_i++)
		{
#pragma HLS PIPELINE
			if(temp_i<pd[rec_a])
			{
				offset=pi[rec_a]+temp_i;
				temp_node=pn[offset];
				if(C1.exist(temp_node,PosC1)||P1.exist(temp_node,PosP1))
				{
					G_index[temp_node]=1;
				}
			}
		}
		if(temp_vector.size)
		{
			PIVOT_3_8:
			for(temp_i=MAX-1;temp_i>=0;temp_i--)
			{
#pragma HLS PIPELINE
				if(temp_i<C2.size)
				{
					temp_node=C2.data[temp_i];
					if(degInP[temp_node]==temp_size2) continue;
					if(degInP[temp_node]<temp_size2-k||degInG[temp_node]<temp_stand1-k)
					{
						C2.remove(temp_node,PosC2);
						PIVOT_3_8_1:
						for(temp_j=0;temp_j<MAX_DEGREE;temp_j++)
						{
							if(temp_j<pd[temp_node])
							{
								offset=pi[temp_node]+temp_j;
								temp_node2=pn[offset];
								degInG[temp_node2]--;
							}
						}
					}
					else
					{
						PIVOT_3_8_2:
						for(nodeT tv=0;tv<MAX;tv++)
						{
							if(tv<temp_vector.size)
							{
								temp_nodei=temp_vector.data[tv];
								if(!bit_G[temp_node][temp_nodei])
								{
									C2.remove(temp_node,PosC2);
									PIVOT_6_2_1:
									for(temp_j=0;temp_j<MAX_DEGREE;temp_j++)
									{
										if(temp_j<pd[temp_node]){
											offset=pi[temp_node]+temp_j;
											temp_node2=pn[offset];
											degInG[temp_node2]--;
										}
									}
									break;
								}
							}
						}
						if(C2.exist(temp_node,PosC2))
						{
							count=0;
							PIVOT_3_8_3:
							for(temp_j=0;temp_j<MAX_DEGREE;temp_j++)
							{
								if(temp_j<pd[temp_node])
								{
									offset=pi[temp_node]+temp_j;
									temp_node2=pn[offset];
									if(G_index[temp_node2]){
										count++;
									}
								}
							}
							if(count<temp_stand1-2*k)
							{
								C2.remove(temp_node,PosC2);
								PIVOT_6_4:
								for(temp_j=0;temp_j<MAX_DEGREE;temp_j++)
								{
									if(temp_j<pd[temp_node])
									{
										offset=pi[temp_node]+temp_j;
										temp_node2=pn[offset];
										degInG[temp_node2]--;
									}
								}
							}
						}
					}
				}
			}
			PIVOT_3_9:
			for(temp_i=MAX-1;temp_i>=0;temp_i--)
			{
#pragma HLS PIPELINE
				if(temp_i<E2.size)
				{
					temp_node=E2.data[temp_i];
					if(degInP[temp_node]==temp_size2) continue;
					if(degInP[temp_node]<temp_size2-k||degInG[temp_node]<temp_stand1-k)
					{
						E2.remove(temp_node,PosE2);
					}
					else
					{
						PIVOT_3_9_1:
						for(nodeT tv=0;tv<MAX;tv++)
						{
							if(tv<temp_vector.size)
							{
								temp_nodei=temp_vector.data[tv];
								if(!bit_G[temp_node][temp_nodei])
								{
									E2.remove(temp_node,PosE2);
									break;
								}
							}
						}
					}
				}
			}
		}
		PIVOT_3_10:
		for(temp_i=0;temp_i<MAX_DEGREE;temp_i++)
		{
#pragma HLS PIPELINE
			if(temp_i<pd[rec_a])
			{
				offset=pi[rec_a]+temp_i;
				temp_node=pn[offset];
				G_index[temp_node]=0;
			}
		}
		if(C2.size+P2.size<temp_stand2)
		{break;}
		PIVOT_3_11:
		while(idx<doing.size&&!C2.exist(doing.data[idx],PosC2))
		{
#pragma HLS PIPELINE
p++;idx++;}
		if(idx==doing.size) break;
		rec_b=doing.data[idx];
		C2.remove(rec_b,PosC2);
		E2.add(rec_b,PosE2);
		PIVOT_3_12:
		for(temp_i=0;temp_i<MAX_DEGREE;temp_i++)
		{
#pragma HLS PIPELINE
			if(temp_i<pd[rec_b])
			{
				offset=pi[rec_b]+temp_i;
				temp_node=pn[offset];
				degInG[temp_node]--;
			}
		}
		branch new_temp;
		if(pivot<pb) new_temp.setvalue(P1,P2,C1,C2,E1,E2);
		else new_temp.setvalue(P2,P1,C2,C1,E2,E1);
		if(verification(new_temp,degInG,left,right,last,mbp,maximum))
		{store_branch(new_temp,ddr,ddrsize,buffer_in_bram,bramsize);}
		E2.remove(rec_b,PosE2);
		rec_a=rec_b;
		idx++;
	}
	if(C1.size+P1.size>=temp_stand1&&C2.size+P2.size>=temp_stand2)
	{
		count=0;
		PIVOT_4:
		for(;idx<doing.size;idx++)
		{
#pragma HLS PIPELINE
			temp_node=doing.data[idx];
			if(C2.exist(temp_node,PosC2)) count++;
		}
		if(C2.size+P2.size-count>=temp_stand2)
		{
			PIVOT_4_1:
			for(;idx<doing.size;idx++)
			{
				temp_node=doing.data[idx];
				if(C2.exist(temp_node,PosC2))
				{
					C2.remove(temp_node,PosC2);
					PIVOT_4_1_1:
					for(temp_i=0;temp_i<MAX_DEGREE;temp_i++)
					{
#pragma HLS PIPELINE
						if(temp_i<pd[temp_node])
						{
							offset=pi[temp_node]+temp_i;
							temp_node2=pn[offset];
							degInG[temp_node2]--;
						}
					}
				}
			}
			if(pivot<pb) new_temp.setvalue(P1,P2,C1,C2,E1,E2);
			else new_temp.setvalue(P2,P1,C2,C1,E2,E1);
			if(verification(new_temp,degInG,left,right,last,mbp,maximum))
			{store_branch(new_temp,ddr,ddrsize,buffer_in_bram,bramsize);}
		}
	}
}
void expand2(branch * ddr,nodeT & ddrsize,branch * buffer_in_bram,nodeT & bramsize,
		degreeT * degInG,degreeT * degInP,nodeT * G_index,
		nodeT * pi,nodeT *pn,degreeT * pd,nodeT pgs,nodeT pb,nodeT left,nodeT right,nodeT k,nodeT last,branch & mbp,nodeT & maximum,
		nodeT pivot,parvector P1,parvector P2,parvector C1,parvector C2,parvector E1,parvector E2,nodeT temp_stand1,nodeT temp_stand2,nodeT * PosP1,nodeT * PosP2,nodeT * PosC1,nodeT * PosC2,nodeT * PosE1,nodeT * PosE2)
{
#pragma HLS INLINE off
	nodeT temp_i;
	nodeT temp_j;
	nodeT temp_k;
	nodeT count;
	nodeT temp_value;
	nodeT temp_size1;
	nodeT temp_size2;
	nodeT temp_node;
	nodeT temp_node2;
	nodeT temp_node3;
	nodeT temp_nodei;
	parvector temp_vector;
	temp_vector.init();
	nodeT offset;
	branch new_temp;
	C1.remove(pivot,PosC1);
	E1.add(pivot,PosE1);
	CPIVOT_1:
	for(temp_i=0;temp_i<MAX_DEGREE;temp_i++)
	{
#pragma HLS PIPELINE
		if(temp_i<pd[pivot])
		{
			offset=pi[pivot]+temp_i;
			temp_node=pn[offset];
			degInG[temp_node]--;
		}
	}
	if(pivot<pb) new_temp.setvalue(P1,P2,C1,C2,E1,E2);
	else new_temp.setvalue(P2,P1,C2,C1,E2,E1);
	if(verification(new_temp,degInG,left,right,last,mbp,maximum))
	{store_branch(new_temp,ddr,ddrsize,buffer_in_bram,bramsize);}
	P1.add(pivot,PosP1);
	CPIVOT_2:
	for(temp_i=0;temp_i<MAX_DEGREE;temp_i++)
	{
#pragma HLS PIPELINE
		if(temp_i<pd[pivot])
		{
			offset=pi[pivot]+temp_i;
			temp_node=pn[offset];
			degInG[temp_node]++;
			degInP[temp_node]++;
		}
	}
	if(degInP[pivot]==P2.size-k)
	{
		CPIVOT_3:
		for(temp_i=MAX-1;temp_i>=0;temp_i--)
		{
#pragma HLS PIPELINE
			if(temp_i<C2.size)
			{
				temp_node=C2.data[temp_i];
				if(!bit_G[pivot][temp_node]||degInG[temp_node]<temp_stand1-k)
				{
					C2.remove(temp_node,PosC2);
					CPIVOT_3_1:
					for(temp_j=0;temp_j<MAX_DEGREE;temp_j++)
					{
						if(temp_j<pd[temp_node])
						{
							offset=pi[temp_node]+temp_j;
							temp_node2=pn[offset];
							degInG[temp_node2]--;
						}
					}
				}
			}
		}
		CPIVOT_4:
		for(temp_i=MAX-1;temp_i>=0;temp_i--)
		{
#pragma HLS PIPELINE
			if(temp_i<E2.size)
			{
				temp_node=E2.data[temp_i];
				if(!bit_G[pivot][temp_node]||degInG[temp_node]<temp_stand1-k)
				{
					E2.remove(temp_node,PosE2);
				}
			}
		}
	}
	else
	{
		temp_size1=P1.size-k;
		CPIVOT_5:
		for(temp_i=MAX-1;temp_i>=0;temp_i--)
		{
#pragma HLS PIPELINE
			if(temp_i<C2.size)
			{
				temp_node=C2.data[temp_i];
				if(degInP[temp_node]<temp_size1||degInG[temp_node]<temp_stand1-k)
				{
					C2.remove(temp_node,PosC2);
					CPIVOT_5_1:
					for(temp_j=0;temp_j<MAX_DEGREE;temp_j++)
					{
						if(temp_j<pd[temp_node])
						{
							offset=pi[temp_node]+temp_j;
							temp_node2=pn[offset];
							degInG[temp_node2]--;
						}
					}
				}
			}
		}
		CPIVOT_6:
		for(temp_i=MAX-1;temp_i>=0;temp_i--)
		{
#pragma HLS PIPELINE
			if(temp_i<E2.size)
			{
				temp_node=E2.data[temp_i];
				if(degInP[temp_node]<temp_size1||degInG[temp_node]<temp_stand1-k)
				{
					E2.remove(temp_node,PosE2);
				}
			}
		}
	}
	CPIVOT_7:
	for(temp_i=0;temp_i<MAX_DEGREE;temp_i++)
	{
#pragma HLS PIPELINE
		if(temp_i<pd[pivot])
		{
			offset=pi[pivot]+temp_i;
			temp_node=pn[offset];
			if(C2.exist(temp_node,PosC2)||P2.exist(temp_node,PosP2))
			{
				G_index[temp_node]=1;
			}
		}
	}
	temp_vector.init();
	temp_size1=P1.size;
	temp_size2=P2.size;
	CPIVOT_8:
	for(temp_i=MAX-1;temp_i>=0;temp_i--)
	{
#pragma HLS PIPELINE
		if(temp_i<P2.size)
		{
			temp_node=P2.data[temp_i];
			if(degInP[temp_node]==temp_size1-k&&!bit_G[pivot][temp_node])
			{
				temp_vector.insert(temp_node);
			}
		}
	}
	if(temp_vector.size&&C2.size+P2.size>=temp_stand2)
	{
		CPIVOT_9:
		for(temp_i=MAX-1;temp_i>=0;temp_i--)
		{
#pragma HLS PIPELINE
			if(temp_i<C1.size)
			{
				temp_node=C1.data[temp_i];
				if(degInP[temp_node]==temp_size2) continue;
				if(degInP[temp_node]<temp_size2-k||degInG[temp_node]<temp_stand2-k)
				{
					C1.remove(temp_node,PosC1);
					CPIVOT_9_1:
					for(temp_j=0;temp_j<MAX_DEGREE;temp_j++)
					{
						if(temp_j<pd[temp_node])
						{
							offset=pi[temp_node]+temp_j;
							temp_node2=pn[offset];
							degInG[temp_node2]--;
						}
					}
				}
				else
				{
					CPIVOT_9_2:
					for(nodeT tv=0;tv<MAX;tv++)
					{
						if(tv<temp_vector.size)
						{
							temp_nodei=temp_vector.data[tv];
							if(!bit_G[temp_node][temp_nodei])
							{
								C1.remove(temp_node,PosC1);
								CPIVOT_9_2_1:
								for(temp_j=0;temp_j<MAX_DEGREE;temp_j++)
								{
									if(temp_j<pd[temp_node])
									{
										offset=pi[temp_node]+temp_j;
										temp_node2=pn[offset];
										degInG[temp_node2]--;
									}
								}
								break;
							}
						}
					}
					if(C1.exist(temp_node,PosC1))
					{
						count=0;
						CPIVOT_9_3:
						for(temp_j=0;temp_j<MAX_DEGREE;temp_j++)
						{
							if(temp_j<pd[temp_node])
							{
								offset=pi[temp_node]+temp_j;
								temp_node2=pn[offset];
								if(G_index[temp_node2]) count++;
							}
						}
						if(count<temp_stand2-2*k)
						{
							C1.remove(temp_node,PosC1);
							CPIVOT_9_4:
							for(temp_j=0;temp_j<MAX_DEGREE;temp_j++)
							{
								if(temp_j<pd[temp_node])
								{
									offset=pi[temp_node]+temp_j;
									temp_node2=pn[offset];
									degInG[temp_node2]--;
								}
							}
						}
					}
				}
			}
		}
		CPIVOT_10:
		for(temp_i=MAX-1;temp_i>=0;temp_i--)
		{
#pragma HLS PIPELINE
			if(temp_i<E1.size)
			{
				temp_node=E1.data[temp_i];
				if(degInP[temp_node]==temp_size2) continue;
				if(degInP[temp_node]<temp_size2-k||degInG[temp_node]<temp_stand2-k)
				{
					E1.remove(temp_node,PosE1);
				}
				else
				{
					CPIVOT_10_1:
					for(nodeT tv=0;tv<MAX;tv++)
					{
						if(tv<temp_vector.size)
						{
							temp_nodei=temp_vector.data[tv];
							if(!bit_G[temp_node][temp_nodei])
							{
								E1.remove(temp_node,PosE1);
								break;
							}
						}
					}
				}
			}
		}
	}
	if(pivot<pb) new_temp.setvalue(P1,P2,C1,C2,E1,E2);
	else new_temp.setvalue(P2,P1,C2,C1,E2,E1);
	if(verification(new_temp,degInG,left,right,last,mbp,maximum))
	{store_branch(new_temp,ddr,ddrsize,buffer_in_bram,bramsize);}
}
void mbs(branch * ddr,nodeT * pi,nodeT *pn,degreeT * pd,nodeT pgs,nodeT pb,nodeT left,nodeT right,nodeT k,nodeT last,branch & mbp,nodeT & maximum)
{
#pragma HLS INLINE off
#ifdef __SYNTHESIS__
	branch batch_for_process[BATCH_SIZE];
	branch buffer_in_bram[BRAM_SIZE];
#else
	branch * batch_for_process=new branch[BATCH_SIZE];
	branch * buffer_in_bram=new branch[BRAM_SIZE];
#endif
	nodeT batchsize=0;
	nodeT bramsize=0;
	nodeT ddrsize=0;
	//define temp
#ifdef __SYNTHESIS__
	degreeT degInP[GRAPH_SIZE];
	degreeT degInG[GRAPH_SIZE];
	nodeT G_index[GRAPH_SIZE];
	nodeT Pos_PL[GRAPH_SIZE];
	nodeT Pos_PR[GRAPH_SIZE];
	nodeT Pos_CL[GRAPH_SIZE];
	nodeT Pos_CR[GRAPH_SIZE];
	nodeT Pos_EL[GRAPH_SIZE];
	nodeT Pos_ER[GRAPH_SIZE];
#else
	degreeT * degInP=new degreeT[GRAPH_SIZE];
	degreeT * degInG=new degreeT[GRAPH_SIZE];
	nodeT * G_index=new nodeT[GRAPH_SIZE];
	nodeT * Pos_PL=new nodeT [GRAPH_SIZE];
	nodeT * Pos_PR=new nodeT [GRAPH_SIZE];
	nodeT * Pos_CL=new nodeT [GRAPH_SIZE];
	nodeT * Pos_CR=new nodeT [GRAPH_SIZE];
	nodeT * Pos_EL=new nodeT [GRAPH_SIZE];
	nodeT * Pos_ER=new nodeT [GRAPH_SIZE];
#endif
#ifdef DEBUG_REDUCTION
	cout<<" in mbs module. "<<endl;
#endif
	nodeT temp_i;
	nodeT temp_j;
	nodeT temp_k;
	nodeT count;
	nodeT pivot;
	nodeT temp_value;
	nodeT temp_size1;
	nodeT temp_size2;
	nodeT temp_node;
	nodeT temp_node2;
	nodeT temp_node3;
	nodeT temp_nodei;
	parvector temp_vector;
	temp_vector.init();
	nodeT offset;
	bool flag=true;
	INIT_BIT_G_I_1:
	for(int init_i=0;init_i<GRAPH_SIZE;init_i++)
	{
#pragma HLS PIPELINE
		INIT_BIT_G_J_1:
		for(int init_j=0;init_j<GRAPH_SIZE;init_j++)
		{
			bit_G[init_i][init_j]=0;
		}
	}
	LOAD_BIT_G_I:
	for(nodeT init_i=0;init_i<pb;init_i++)
	{
#pragma HLS PIPELINE
		LOAD_BIT_G_J:
		for(nodeT init_j=0;init_j<MAX_DEGREE;init_j++)
		{
			if(init_j<pd[init_i])
			{
				offset=pi[init_i]+init_j;
				temp_node=pn[offset];
				bit_G[init_i][temp_node]=true;
				bit_G[temp_node][init_i]=true;
			}
		}
	}
	EACH_NODE:
	for(nodeT i=0;i<pb;i++)
	{
#ifdef DEBUG_REDUCTION
		cout<<" before init and deal node:  "<<i<<endl;
#endif
		init_deal_space(degInP,degInG,G_index,Pos_PL,Pos_PR,Pos_CL,Pos_CR,Pos_EL,Pos_ER);
		//init the first branch
		branch b0;
		b0.init();
		b0.P_L.add(i,Pos_PL);
#ifdef DEBUG_REDUCTION
		cout<<" after init and then go to MBS init 1 loop "<<endl;
#endif
		MBS_INIT_1_LOOP:
		for(temp_i=0;temp_i<MAX_DEGREE;temp_i++)
		{
#pragma HLS PIPELINE
			if(temp_i<pd[i])
			{
				offset=pi[i]+temp_i;
				temp_node=pn[offset];
				count=0;
				MBS_INIT_1_1_LOOP:
				for(temp_j=0;temp_j<MAX_DEGREE;temp_j++)
				{
					if(temp_j<pd[temp_node])
					{
						offset=pi[temp_node]+temp_j;
						temp_node2=pn[offset];
						if(temp_node2>=i) count++;
					}
				}
				if(count>=left-k)
				{
					degInG[temp_node]++;
					degInP[temp_node]++;
					b0.C_R.add(temp_node,Pos_CR);
					MBS_INIT_1_2_LOOP:
					for(temp_j=0;temp_j<MAX_DEGREE;temp_j++)
					{
						if(temp_j<pd[temp_node])
						{
							offset=pi[temp_node]+temp_j;
							temp_node2=pn[offset];
							degInG[temp_node2]++;
						}
					}
				}
			}
		}
#ifdef DEBUG_REDUCTION
		cout<<" go to MBS init 2 loop "<<endl;
#endif
		MBS_INIT_2_LOOP:
		for(temp_i=0;temp_i<MAX_DEGREE;temp_i++)
		{
			if(temp_i<pd[i])
			{
				offset=pi[i]+temp_i;
				temp_node=pn[offset];
				MBS_INIT_2_1_LOOP:
				for(temp_j=0;temp_j<MAX_DEGREE;temp_j++)
				{
#pragma HLS PIPELINE
					if(temp_j<pd[temp_node])
					{
						offset=pi[temp_node]+temp_j;
						temp_node2=pn[offset];

						if(temp_node2>i&&degInG[temp_node2]>=right-2*k&&!b0.C_L.exist(temp_node2,Pos_CL))
						{
							b0.C_L.add(temp_node2,Pos_CL);
							MBS_INIT_2_1_1_LOOP:
							for(temp_k=0;temp_k<MAX_DEGREE;temp_k++)
							{
								if(temp_k<pd[temp_node2])
								{
									offset=pi[temp_node2]+temp_k;
									temp_node3=pn[offset];
									degInG[temp_node3]++;
								}
							}
						}
						if(temp_node2<i&&degInG[temp_node2]>=right-2*k&&!b0.E_L.exist(temp_node2,Pos_EL))
						{
							b0.E_L.add(temp_node2,Pos_EL);
						}
					}
				}
			}
		}
#ifdef DEBUG_REDUCTION
		cout<<" go to MBS init prune 1"<<endl;
#endif
		MBS_INIT_PRUNE_1:
		for(temp_i=0;temp_i<GRAPH_SIZE;temp_i++)
		{
			if(temp_i>=pb&&temp_i<pgs)
			{
				bool temp_exist=b0.C_R.exist(temp_i,Pos_CR);
				if(degInG[temp_i]>=left-k&&!temp_exist)
				{
					b0.C_R.add(temp_i,Pos_CR);
					MBS_INIT_PRUNE_1_1:
					for(temp_j=0;temp_j<MAX_DEGREE;temp_j++)
					{
#pragma HLS PIPELINE
						if(temp_j<pd[temp_i])
						{
							offset=pi[temp_i]+temp_j;
							temp_node=pn[offset];
							degInG[temp_node]++;
						}
					}
				}
				if(degInG[temp_i]<left-k&&temp_exist)
				{
					b0.C_R.remove(temp_i,Pos_CR);
					MBS_INIT_PRUNE_1_2:
					for(temp_j=0;temp_j<MAX_DEGREE;temp_j++)
					{
#pragma HLS PIPELINE
						if(temp_j<pd[temp_i])
						{
							offset=pi[temp_i]+temp_j;
							temp_node=pn[offset];
							degInG[temp_node]--;
							if(degInG[temp_node]<right-k&&b0.C_L.exist(temp_node,Pos_CL))
							{
								b0.C_L.remove(temp_node,Pos_CL);
								MBS_INIT_PRUNE_1_2_1:
								for(temp_k=0;temp_k<MAX_DEGREE;temp_k++)
								{
									if(temp_k<pd[temp_node])
									{
										offset=pi[temp_node]+temp_k;
										temp_node2=pn[offset];
										degInG[temp_node2]--;
									}
								}
							}
						}
					}
				}
			}
		}
#ifdef DEBUG_REDUCTION
		cout<<" go to MBS init prune 2 "<<endl;
#endif
		MBS_INIT_PRUNE_2:
		for(temp_i=0;temp_i<GRAPH_SIZE;temp_i++)
		{
			if(temp_i>=pb&&temp_i<pgs)
			{
				bool temp_exist=b0.C_R.exist(temp_i,Pos_CR);
				if(degInG[temp_i]>=left-k&&!temp_exist)
				{
					b0.C_R.add(temp_i,Pos_CR);
					MBS_INIT_PRUNE_2_1:
					for(temp_j=0;temp_j<MAX_DEGREE;temp_j++)
					{
#pragma HLS PIPELINE
						if(temp_j<pd[temp_i])
						{
							offset=pi[temp_i]+temp_j;
							temp_node=pn[offset];
							degInG[temp_node]++;
						}
					}
				}
				if(degInG[temp_i]<left-k&&temp_exist)
				{
					b0.C_R.remove(temp_i,Pos_CR);
					MBS_INIT_PRUNE_2_2:
					for(temp_j=0;temp_j<MAX_DEGREE;temp_j++)
					{
#pragma HLS PIPELINE
						if(temp_j<pd[temp_i])
						{
							offset=pi[temp_i]+temp_j;
							temp_node=pn[offset];
							degInG[temp_node]--;
							if(degInG[temp_node]<right-k&&b0.C_L.exist(temp_node,Pos_CL))
							{
								b0.C_L.remove(temp_node,Pos_CL);
								MBS_INIT_PRUNE_2_2_1:
								for(temp_k=0;temp_k<MAX_DEGREE;temp_k++)
								{
									if(temp_k<pd[temp_node])
									{
										offset=pi[temp_node]+temp_k;
										temp_node2=pn[offset];
										degInG[temp_node2]--;
									}
								}
							}
						}
					}
				}
			}
		}
#ifdef DEBUG_REDUCTION
		cout<<" go to MBS init prune 3"<<endl;
#endif
		MBS_INIT_PRUNE_3:
		for(temp_i=0;temp_i<GRAPH_SIZE;temp_i++)
		{
			if(temp_i>=pb&&temp_i<pgs)
			{
				bool temp_exist=b0.C_R.exist(temp_i,Pos_CR);
				if(degInG[temp_i]>=left-k&&!temp_exist)
				{
					b0.C_R.add(temp_i,Pos_CR);
					MBS_INIT_PRUNE_3_1:
					for(temp_j=0;temp_j<MAX_DEGREE;temp_j++)
					{
#pragma HLS PIPELINE
						if(temp_j<pd[temp_i])
						{
							offset=pi[temp_i]+temp_j;
							temp_node=pn[offset];
							degInG[temp_node]++;
						}
					}
				}
				if(degInG[temp_i]<left-k&&temp_exist)
				{
					b0.C_R.remove(temp_i,Pos_CR);
					MBS_INIT_PRUNE_3_2:
					for(temp_j=0;temp_j<MAX_DEGREE;temp_j++)
					{
#pragma HLS PIPELINE
						if(temp_j<pd[temp_i])
						{
							offset=pi[temp_i]+temp_j;
							temp_node=pn[offset];
							degInG[temp_node]--;
							if(degInG[temp_node]<right-k&&b0.C_L.exist(temp_node,Pos_CL))
							{
								b0.C_L.remove(temp_node,Pos_CL);
								MBS_INIT_PRUNE_3_2_1:
								for(temp_k=0;temp_k<MAX_DEGREE;temp_k++)
								{
									if(temp_k<pd[temp_node])
									{
										offset=pi[temp_node]+temp_k;
										temp_node2=pn[offset];
										degInG[temp_node2]--;
									}
								}
							}
						}
					}
				}
			}
		}

		flag=true;
		FLAG_JUDGE:
		for(temp_i=MAX-1;temp_i>=0;temp_i--)
		{
			if(temp_i<b0.E_L.size&&flag)
			{
				temp_node=b0.E_L.data[temp_i];
				if(degInG[temp_node]==b0.C_R.size)
				{
					flag=false;
				}
			}
		}
		batchsize=0;
		bramsize=0;
		ddrsize=0;
		batch_for_process[batchsize]=b0;
		batchsize++;
		flag=verification(b0,degInG,left,right,last,mbp,maximum);


		EXPANSION_AND_VERIFICATION:
		while(flag&&batchsize)
		{

			BATCH_DEAL:
			for(nodeT temp_b=0;temp_b<BATCH_SIZE;temp_b++)
			{
				if(temp_b<batchsize)
				{
					bool judge_return=false;
					branch new_temp;
					init_deal_space(degInP,degInG,G_index,Pos_PL,Pos_PR,Pos_CL,Pos_CR,Pos_EL,Pos_ER);
					branch temp_branch=batch_for_process[temp_b];
					CONSTRUCT_DEGREE_PL:
					for(temp_j=0;temp_j<MAX;temp_j++)
					{
#pragma HLS PIPELINE
						if(temp_j<temp_branch.P_L.size)
						{
							temp_node=temp_branch.P_L.data[temp_j];
							CONSTRUCT_DEGREE_PL1:
							for(temp_k=0;temp_k<MAX_DEGREE;temp_k++)
							{
								if(temp_k<pd[temp_node])
								{
									offset=pi[temp_node]+temp_k;
									temp_node2=pn[offset];
									degInP[temp_node2]++;
									degInG[temp_node2]++;
								}
							}
						}
					}
					CONSTRUCT_DEGREE_PR:
					for(temp_j=0;temp_j<MAX;temp_j++)
					{
#pragma HLS PIPELINE
						if(temp_j<temp_branch.P_R.size)
						{
							temp_node=temp_branch.P_R.data[temp_j];
							CONSTRUCT_DEGREE_PR2:
							for(temp_k=0;temp_k<MAX_DEGREE;temp_k++)
							{
								if(temp_k<pd[temp_node])
								{
									offset=pi[temp_node]+temp_k;
									temp_node2=pn[offset];
									degInP[temp_node2]++;
									degInG[temp_node2]++;
								}
							}
						}
					}
					CONSTRUCT_DEGREE_CL:
					for(temp_j=0;temp_j<MAX;temp_j++)
					{
#pragma HLS PIPELINE
						if(temp_j<temp_branch.C_L.size)
						{
							temp_node=temp_branch.C_L.data[temp_j];
							CONSTRUCT_DEGREE_CL2:
							for(temp_k=0;temp_k<MAX_DEGREE;temp_k++)
							{
								if(temp_k<pd[temp_node])
								{
									offset=pi[temp_node]+temp_k;
									temp_node2=pn[offset];
									degInG[temp_node2]++;
								}
							}
						}
					}
					CONSTRUCT_DEGREE_CR:
					for(temp_j=0;temp_j<MAX;temp_j++)
					{
#pragma HLS PIPELINE
						if(temp_j<temp_branch.C_R.size)
						{
							temp_node=temp_branch.C_R.data[temp_j];
							CONSTRUCT_DEGREE_CR2:
							for(temp_k=0;temp_k<MAX_DEGREE;temp_k++)
							{
								if(temp_k<pd[temp_node])
								{
									offset=pi[temp_node]+temp_k;
									temp_node2=pn[offset];
									degInG[temp_node2]++;
								}
							}
						}
					}
					loadpos(temp_branch.P_L,Pos_PL);
					loadpos(temp_branch.P_R,Pos_PR);
					loadpos(temp_branch.C_L,Pos_CL);
					loadpos(temp_branch.C_R,Pos_CR);
					loadpos(temp_branch.E_L,Pos_EL);
					loadpos(temp_branch.E_R,Pos_ER);
					pivot=-1;
					temp_value=0;
					temp_node2=pgs;
					temp_j=temp_branch.P_L.size+temp_branch.C_L.size;
					temp_k=temp_branch.P_R.size+temp_branch.C_R.size;
					FIND_PIVOT_PL:
					for(temp_i=MAX-1;temp_i>=0;temp_i--)
					{
#pragma HLS PIPELINE
						if(temp_i<temp_branch.P_L.size)
						{
							temp_node=temp_branch.P_L.data[temp_i];
							if(temp_k-degInG[temp_node]>temp_value)
							{
								temp_value=temp_k-degInG[temp_node];
								pivot=temp_node;
							}
							if(temp_node2>degInG[temp_node])
							{
								temp_node2=degInG[temp_node];
							}
						}
					}
					if(temp_node2+k<right) continue;
					temp_size1=pgs;
					FIND_PIVOT_PR:
					for(temp_i=MAX-1;temp_i>=0;temp_i--)
					{
#pragma HLS PIPELINE
						if(temp_i<temp_branch.P_R.size)
						{
							temp_node=temp_branch.P_R.data[temp_i];
							if(temp_j-degInG[temp_node]>temp_value)
							{
								temp_value=temp_j-degInG[temp_node];
								pivot=temp_node;
							}
							if(temp_size1>degInG[temp_node])
							{
								temp_size1=degInG[temp_node];
							}
						}
					}
					if(temp_size1+k<left||(temp_node2+k)*(temp_size1+k)<=maximum+2*k) continue;
					if(pivot>=0&&temp_value>k)
					{
						if(pivot<pb)
						{
							expand1(ddr,ddrsize,buffer_in_bram,bramsize,degInG,degInP,G_index,pi,pn,pd,pgs,pb,left,right,k,last,mbp,maximum,pivot,
									temp_branch.P_L,temp_branch.P_R,temp_branch.C_L,temp_branch.C_R,temp_branch.E_L,temp_branch.E_R,left,right,
									Pos_PL,Pos_PR,Pos_CL,Pos_CR,Pos_EL,Pos_ER);
						}
						else
						{
							expand1(ddr,ddrsize,buffer_in_bram,bramsize,degInG,degInP,G_index,pi,pn,pd,pgs,pb,left,right,k,last,mbp,maximum,pivot,
									temp_branch.P_R,temp_branch.P_L,temp_branch.C_R,temp_branch.C_L,temp_branch.E_R,temp_branch.E_L,right,left,
									Pos_PR,Pos_PL,Pos_CR,Pos_CL,Pos_ER,Pos_EL);
						}

					}
					else
					{//not find pivot in P then find pivot in C
						FIND_PIVOT_CL:
						for(temp_i=MAX-1;temp_i>=0;temp_i--)
						{
#pragma HLS PIPELINE
							if(temp_i<temp_branch.C_L.size)
							{
								temp_node=temp_branch.C_L.data[temp_i];
								if(temp_k-degInG[temp_node]>temp_value)
								{
									temp_value=temp_k-degInG[temp_node];
									pivot=temp_node;
								}
							}
						}
						FIND_PIVOT_CR:
						for(temp_i=MAX-1;temp_i>=0;temp_i--)
						{
#pragma HLS PIPELINE
							if(temp_i<temp_branch.C_R.size)
							{
								temp_node=temp_branch.C_R.data[temp_i];
								if(temp_j-degInG[temp_node]>temp_value)
								{
									temp_value=temp_j-degInG[temp_node];
									pivot=temp_node;
								}
							}
						}
						//if not has a pivot
						if(temp_value<=k)
						{
							temp_value=0;
							NO_PIVOT_1:
							for(temp_i=MAX-1;temp_i>=0;temp_i--)
							{
#pragma HLS PIPELINE
								if(temp_i<temp_branch.P_L.size)
								{
									temp_node=temp_branch.P_L.data[temp_i];
									temp_value+=degInG[temp_node];
								}
							}
							NO_PIVOT_2:
							for(temp_i=MAX-1;temp_i>=0;temp_i--)
							{
#pragma HLS PIPELINE
								if(temp_i<temp_branch.C_L.size)
								{
									temp_node=temp_branch.C_L.data[temp_i];
									temp_value+=degInG[temp_node];
								}
							}
							if(maximum<temp_value)
							{
								maximum=temp_value;
								mbp=temp_branch;
							}
							continue;
						}
						if(pivot>=0&&temp_value>k)
						{
							if(pivot<pb)
							{
								expand2(ddr,ddrsize,buffer_in_bram,bramsize,degInG,degInP,G_index,pi,pn,pd,pgs,pb,left,right,k,last,mbp,maximum,pivot,
										temp_branch.P_L,temp_branch.P_R,temp_branch.C_L,temp_branch.C_R,temp_branch.E_L,temp_branch.E_R,left,right,
										Pos_PL,Pos_PR,Pos_CL,Pos_CR,Pos_EL,Pos_ER);
							}
							else
							{
								expand2(ddr,ddrsize,buffer_in_bram,bramsize,degInG,degInP,G_index,pi,pn,pd,pgs,pb,left,right,k,last,mbp,maximum,pivot,
										temp_branch.P_R,temp_branch.P_L,temp_branch.C_R,temp_branch.C_L,temp_branch.E_R,temp_branch.E_L,right,left,
										Pos_PR,Pos_PL,Pos_CR,Pos_CL,Pos_ER,Pos_EL);
							}
						}
					}

				}
			}
			get_branch(batch_for_process,batchsize,buffer_in_bram,bramsize,ddr,ddrsize);
		}
	}
}
void record(branch & b,nodeT * g)
{
	branch temp;
	nodeT index;
	nodeT temp_num;
	temp.init();
	RECORD_1:
	for(nodeT i=0;i<MAX;i++)
	{
#pragma HLS PIPELINE
		if(i<b.P_L.size)
		{
			index=b.P_L.data[i];
			temp_num=g[index];
			temp.P_L.insert(temp_num);
		}
	}
	RECORD_2:
	for(nodeT i=0;i<MAX;i++)
	{
#pragma HLS PIPELINE
		if(i<b.C_L.size)
		{
			index=b.C_L.data[i];
			temp_num=g[index];
			temp.C_L.insert(temp_num);
		}
	}
	RECORD_3:
	for(nodeT i=0;i<MAX;i++)
	{
#pragma HLS PIPELINE
		if(i<b.P_R.size)
		{
			index=b.P_R.data[i];
			temp_num=g[index];
			temp.P_R.insert(temp_num);
		}
	}
	RECORD_4:
	for(nodeT i=0;i<MAX;i++)
	{
#pragma HLS PIPELINE
		if(i<b.C_R.size)
		{
			index=b.C_R.data[i];
			temp_num=g[index];
			temp.C_R.insert(temp_num);
		}
	}
	b=temp;
}
extern "C"
{
void Framework(branch * ddr,nodeT * Index,nodeT * Neighbor,degreeT * Degree,nodeT Graph_size,nodeT Bipartite,nodeT K,nodeT Alpha,nodeT Beta,branch * mbp,nodeT * maximum)
{
#pragma HLS INTERFACE m_axi port=ddr offset=slave bundle=gmem1 depth=100
#pragma HLS INTERFACE m_axi port=Index offset=slave bundle=gmem2 depth=45
#pragma HLS INTERFACE m_axi port=Neighbor offset=slave bundle=gmem2 depth=198
#pragma HLS INTERFACE m_axi port=Degree offset=slave bundle=gmem2 depth=44
#pragma HLS INTERFACE m_axi port=mbp offset=slave bundle=gmem2 depth=1
#pragma HLS INTERFACE m_axi port=maximum offset=slave bundle=gmem2 depth=1

#pragma HLS INTERFACE s_axilite port=Graph_size bundle=control
#pragma HLS INTERFACE s_axilite port=Bipartite bundle=control
#pragma HLS INTERFACE s_axilite port=K bundle=control
#pragma HLS INTERFACE s_axilite port=Alpha bundle=control
#pragma HLS INTERFACE s_axilite port=Beta bundle=control
#pragma HLS INTERFACE s_axilite port=ddr bundle=control
#pragma HLS INTERFACE s_axilite port=Index bundle=control
#pragma HLS INTERFACE s_axilite port=Neighbor bundle=control
#pragma HLS INTERFACE s_axilite port=Degree bundle=control
#pragma HLS INTERFACE s_axilite port=mbp bundle=control
#pragma HLS INTERFACE s_axilite port=maximum bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control
	//parameter load
	nodeT k=K;
	nodeT graph_size=Graph_size;
	nodeT theta_l=Alpha;
	nodeT theta_r=Beta;
	nodeT bipartite=Bipartite;
#ifdef __SYNTHESIS__
	nodeT index[GRAPH_SIZE+1];
	nodeT neighbor[EDGE_SIZE*2];
	degreeT degree[GRAPH_SIZE];
#else
	nodeT * index=new nodeT[GRAPH_SIZE+1];
	nodeT * neighbor=new nodeT[EDGE_SIZE*2];
	degreeT * degree=new degreeT[GRAPH_SIZE];
#endif
	//variable declare and init
	nodeT low_l;
	nodeT low_r;
	nodeT last;
	nodeT maxLdeg=0;
	nodeT maxRdeg=0;
	nodeT current=0;
	nodeT last_current;
	branch maxbp;
	LOAD_INDEX:
	for(nodeT i=0;i<GRAPH_SIZE+1;i++)
	{
#pragma HLS PIPELINE
		index[i]=Index[i];
	}
	LOAD_NEIGHBOR:
	for(nodeT i=0;i<EDGE_SIZE*2;i++)
	{
#pragma HLS PIPELINE
		neighbor[i]=Neighbor[i];
	}
	LOAD_DEGREE:
	for(nodeT i=0;i<GRAPH_SIZE;i++)
	{
#pragma HLS PIPELINE
		degree[i]=Degree[i];
	}
	FIND_DEGREE:
	for(nodeT i=0;i<GRAPH_SIZE;i++)
	{
#pragma HLS PIPELINE
		if(i<bipartite)
		{if(maxLdeg<degree[i]){maxLdeg=degree[i];}}
		else if(i>=bipartite&&i<graph_size)
		{if(maxRdeg<degree[i]){maxRdeg=degree[i];}}
	}
	low_l=maxLdeg;
	PBIE:
	while(low_l>theta_l)
	{
		last=low_l;
		last_current=current;
		if(current/low_l<theta_l)
		{low_r=theta_l;}
		else
		{low_r=current/low_l;}
		if(low_l/2<theta_r)
		{low_l=theta_r;}
		else
		{low_l=low_l/2;}
#ifdef PRINT
		cout<<"now low bound of L is "<<low_l<<" and low bound of R is "<<low_r<<endl;
#endif
		nodeT pgs;
		nodeT pb;
#ifdef __SYNTHESIS__
		nodeT pi[GRAPH_SIZE];
		nodeT pg[GRAPH_SIZE];
		nodeT pn[EDGE_SIZE*2];
		degreeT pd[GRAPH_SIZE];
#else
		nodeT * pi=new nodeT[GRAPH_SIZE];
		nodeT * pg=new nodeT[GRAPH_SIZE];
		nodeT * pn=new nodeT[EDGE_SIZE*2];
		degreeT * pd=new degreeT[GRAPH_SIZE];
#endif
		coredecompositon(index,neighbor,degree,graph_size,bipartite,pg,pi,pn,pd,pgs,pb,low_r-k,low_l-k);
		if(pgs)
		{
			mbs(ddr,pi,pn,pd,pgs,pb,low_l,low_r,k,last,maxbp,current);
			if(current!=last_current)
			{
				record(maxbp,pg);
			}
		}
	}
	mbp[0]=maxbp;
	maximum[0]=current;
}
}
