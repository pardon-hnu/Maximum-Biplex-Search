#include<iostream>
#include<list>
#include<map>
#include<algorithm>
#include <vector>
#include <cmath>
using namespace std;

struct xypair
{
    int low;
    int upper;
    xypair()
    {

    }
    xypair(int l,int u)
    {
        this->low=l;
        this->upper=u;
    }
};
class XYPair
{
    private:
        map<int,xypair> xpair;
        map<int,xypair> ypair;
    public:
        XYPair();
        XYPair(map<int,xypair> p1,map<int,xypair> p2);
        ~XYPair();
        xypair SearchX(int x);
        xypair SearchY(int y);
        int GetMaxX();
        int GetMaxY();
        int GetMinX();
        int GetMinY();
        int getUpperBound(int x,int y);
        int getLowBound(int x,int y,int k);
        void update(int maximum);
        void removeX(int x);
        void removeY(int y);
        void printXYPair();
};
XYPair::XYPair()
{

}
XYPair::XYPair(map<int,xypair> p1,map<int,xypair> p2)
{
    this->xpair=p1;
    this->ypair=p2;
}
XYPair::~XYPair()
{

}
xypair XYPair::SearchX(int x)
{
    return xpair.find(x)->second;
}
xypair XYPair::SearchY(int y)
{
    return ypair.find(y)->second;
}
int XYPair::GetMaxX()
{
    return xpair.rbegin()->first;
}
int XYPair::GetMaxY()
{
    return ypair.rbegin()->first;
}
int XYPair::GetMinX()
{
    return xpair.begin()->first;
}
int XYPair::GetMinY()
{
    return ypair.begin()->first;
}
int XYPair::getUpperBound(int x,int y)
{
return x*y;
}
int XYPair::getLowBound(int x,int y,int k)
{
return (x>y)?(x*y-y*k):(x*y-x*k);
}
void XYPair::update(int maximum)
{
    // cout<<"now remove all pairs that x*y <="<<maximum<<endl;
    map<int,xypair>::iterator it;
    for(it=xpair.begin();it!=xpair.end();)
    {
        int x=it->first;
        xypair p=it->second;
        int new_low=ceil((double)maximum/x);
        // cout<<"now x: "<<x<<" and new low is: "<<new_low<<" old pair: "<<p.low<<"->"<<p.upper<<endl;
        if(new_low>p.upper)
        {
            xpair.erase(it++);
        }
        else
        {  
            it->second.low=new_low;
            it++;
        }
    }
    for(it=ypair.begin();it!=ypair.end();)
    {
        int y=it->first;
        xypair p=it->second;
        int new_low=ceil((double)maximum/y);
        // cout<<"now y: "<<y<<" and new low is: "<<new_low<<" old pair: "<<p.low<<"->"<<p.upper<<endl;
        if(new_low>p.upper)
        {
            ypair.erase(it++);
        }
        else
        {  
            it->second.low=new_low;
            it++;
        }
    }
}
void XYPair::removeX(int x)
{
    map<int,xypair>::iterator it;
    for(it=xpair.begin();it!=xpair.end();)
    {
        if(it->first>=x)
        {
            xpair.erase(it,xpair.end());
            break;
        }
        else
        {
            it++;
        }
    }
    for(it=ypair.begin();it!=ypair.end();)
    {
        if(it->second.low>=x)
        {
            ypair.erase(it++);
        }
        else
        {
            if(it->second.upper>=x)
            {
                it->second.upper=x-1;
            }
            it++;
        }
    }
}
void XYPair::removeY(int y)
{
    // cout<<"remove all pairs that y >= "<<y<<endl;
    // printXYPair();
    map<int,xypair>::iterator it;
    for(it=ypair.begin();it!=ypair.end();)
    {
        if(it->first>=y)
        {
            ypair.erase(it,ypair.end());
            break;
        }
        else
        {
            it++;
        }
    }
    // printXYPair();
    for(it=xpair.begin();it!=xpair.end();)
    {
        if(it->second.low>=y)
        {
            xpair.erase(it++);
        }
        else
        {
            if(it->second.upper>=y)
            {
                it->second.upper=y-1;
            }
            it++;
        }
    }
}
void XYPair::printXYPair()
{
    cout<<"now print x direction: "<<endl;
    map<int,xypair>::iterator it;
    for(it=xpair.begin();it!=xpair.end();it++)
    {
        int x=it->first;
        xypair p=it->second;
        cout<<x<<":"<<p.low<<"-"<<p.upper<<endl;
    }
    cout<<"now print y direction: "<<endl;
    for(it=ypair.begin();it!=ypair.end();it++)
    {
        int y=it->first;
        xypair p=it->second;
        cout<<y<<":"<<p.low<<"-"<<p.upper<<endl;
    }
}