#include <bits/stdc++.h>

using namespace std;

int n, r;


void outputGraph(vector< vector<long long> > &graph)
{
    /*cout << "Outputting a graph" << endl;
    cout << "The graph has " << graph.size() << " vertices" << endl;
    
    for(int i=0; i<graph.size(); i++)
    {
        cout << "Vertex " << i << " has the following neighbors:" << endl;
        for(int neigh : graph[i])
        {
            cout << neigh << " ";
        }
        cout << endl;
    }*/

    vector<int> empVec;
    vector< vector<int> > Levels(r+1,empVec);

    vector<int> dist(graph.size(),1e9);
    dist[0]=0;
    queue<int> q;
    q.push(0);
    Levels[0].push_back(0);
    while(!q.empty())
    {
        int now=q.front();
        q.pop();
        for(int neigh : graph[now])
        {
            if(dist[neigh]<1e9) continue;
            dist[neigh]=dist[now]+1;
            Levels[dist[neigh]].push_back(neigh);
            q.push(neigh);
        }
    }

    vector<int> leaveCounts;

    set<int> empSet;
    vector< set<int> > labelsPerVertex(graph.size(),empSet);

    //queue<int> q;
    int specialCounter=0;
    cout << "Leave counts: ";
    for(int u : Levels[r-2])
    {
        leaveCounts.push_back(graph[u].size()-1);
        cout << graph[u].size()-1 << " ";
        labelsPerVertex[u].insert(specialCounter);
        specialCounter++;
        q.push(u);
    }
    cout << endl;

    while(!q.empty())
    {
        int now=q.front();
        q.pop();
        for(int neigh : graph[now])
        {
            if(dist[neigh]>=dist[now]) continue;
            for(int labelNow : labelsPerVertex[now])
            {
                labelsPerVertex[neigh].insert(labelNow);
            }
            q.push(neigh);
        }
    }


    for(int i=r-3; i>=0; i--)
    {
        for(int u : Levels[i])
        {
            if(graph[u].size()>=3)
            {
                cout << "New vertex with labels ";
                for(int label : labelsPerVertex[u]) cout << label << " ";
                cout << "at distance " << r-1-i << " from the leaves" << endl;
                cout << endl;
            }
        }
    }
}

// command line arguments: n and r
int main(int argc, char ** argv) {

    if(argc!=3)
    {
        fprintf(stderr,"Wrong number of command line arguments!\nExpected n and r\n");
    }
    n=atoi(argv[1]);
    r=atoi(argv[2]);
    ios::sync_with_stdio(false);
    cin.tie(0);
    
    long long inf=1e16;
    vector< vector<long long> > dW(n+2, vector<long long>(r+2,-inf));

    pair<int, int> empPair=make_pair(-1,-1); // capture (t1,t2) or (t-1,0) if corresponding vertex has degree 2
    vector< vector< pair<int, int> > > argMaxdW(n+2, vector< pair<int, int> >(r+2,empPair));
    
    // i=1
    for(int t=1; t<=n; t++)
    {
        dW[t][1]=t*(n-1);
        // argMaxdW is not filled in for this special step
    }

    // i between 2 and r-1
    for(int i=2; i<=r-1; i++)
    {
        for(int t=i; t<=n-r+i-1; t++)
        {
            // option 1
            long long option1=dW[t-1][i-1]+t*(n-t);
            
            // option 2
            long long option2=-inf;
            int bestT1=-1;
            int bestT2=-1;
            for(int t1=1; t1+1+2<=t; t1++)
            {
                int t2=t-2-t1;
                long long now=dW[t1][i-1]+dW[t2][i-1]+(t1+1)*(n-t1-1)+(t2+1)*(n-t2-1);
                if(now>option2)
                {
                    option2=now;
                    bestT1=t1;
                    bestT2=t2;
                }
            }

            if(option1>=option2)
            {
                dW[t][i]=option1;
                argMaxdW[t][i]=make_pair(t-1,0);
            }
            else
            {
                dW[t][i]=option2;
                argMaxdW[t][i]=make_pair(bestT1,bestT2);
            }
        }
    }

    // i=r
    for(int t=r; t<=n-1; t++)
    {
        dW[t][r]=dW[t-1][r-1]+t*(n-t);
        // argMaxdW is not filled in for this special step
    }


    vector<long long> profits;
    vector<long long> weights;

    vector<long long> multipliers;
    vector<int> tis;
    for(int ti=1; ti<=n; ti++)
    {
        for(long long multiplier=1; multiplier*ti<=n-1; multiplier*=2)
        {
            if(dW[ti][r]>=0)
            {
                weights.push_back(multiplier*ti);
                profits.push_back(multiplier*dW[ti][r]);

                multipliers.push_back(multiplier);
                tis.push_back(ti);
            }
        }
    }
    int nbItems=profits.size();
    int capacity=n-1;

    vector< vector<long long> > bestProfit(nbItems+1, vector<long long>(capacity+1,-inf));
    vector< vector< pair<int, int> > > prevState(nbItems+1, vector< pair<int, int> >(capacity+1,empPair));
    
    // base case: no items available
    bestProfit[0][0]=0;

    // recursive case
    for(int i=1; i<=nbItems; i++)
    {
        for(int totWeight=0; totWeight<=capacity; totWeight++)
        {
            // option 1: do not take current item
            long long option1=bestProfit[i-1][totWeight];

            // option 2: take current item
            long long option2=-inf;
            long long prevWeight=totWeight-weights[i-1];
            if(prevWeight>=0)
            {
                option2=profits[i-1]+bestProfit[i-1][prevWeight];
            }

            if(option1>=option2)
            {
                bestProfit[i][totWeight]=option1;
                prevState[i][totWeight]=make_pair(i-1,totWeight);
            }
            else
            {
                bestProfit[i][totWeight]=option2;
                prevState[i][totWeight]=make_pair(i-1,prevWeight);
            }
        }
    }


    // output solution for single centre

    int idx1=nbItems;
    int idx2=n-1;
    long long value=bestProfit[idx1][idx2];

    cout << "Start output for case where the tree has a single centre:" << endl;
    cout << "Maximum Wiener index of a graph with " << n << " vertices and radius " << r << " is: " << value << endl;

    map<int, int> tiAndTimes;
    while(idx1>0)
    {
        pair<int, int> prev=prevState[idx1][idx2];
        int newIdx1=prev.first;
        int newIdx2=prev.second;
        if(newIdx2<idx2)
        {
            int ti=tis[idx1-1];
            tiAndTimes[ti]+=multipliers[idx1-1];
        }
        idx1=newIdx1;
        idx2=newIdx2;
    }

    // output
    for (const auto& pa2 : tiAndTimes) 
    {
        int ti=pa2.first;
        int mult=pa2.second;

        cout << "Take dW(" << n << "," << ti << "," << r << ") " << mult << " times" << endl;


        // build tree

        vector< vector<long long> > tree;
        vector<long long> empVec;
        tree.push_back(empVec);

        // special case i=r: check later how to proceed with this special case
        // TODO

        queue< pair<int, pair<int, int> > > q;
        q.push(make_pair(0,make_pair(ti-1,r-1)));
        while(q.size()>0)
        {
            pair<int, pair<int, int> > pa=q.front();
            q.pop();
            int u=pa.first;

            int t=pa.second.first;
            int i=pa.second.second;

            if(i==1)
            {
                for(int times=1; times<=t; times++)
                {
                    tree.push_back(empVec);
                    int newNode=tree.size()-1;
                    tree[u].push_back(newNode);
                    tree[newNode].push_back(u);
                }
            }
            else // i should be between 2 and r-1
            {
                if(!(2<=i && i<=r-1))
                {
                    cout << "Logical error. The index i should be between 2 and " << r-1 << ", but it was " << i << endl;
                    exit(0);
                }
                pair<int, int> prev=argMaxdW[t][i];
            
                if(prev.second==0) // add 1 neighbor
                {
                    tree.push_back(empVec);
                    int newNode=tree.size()-1;
                    tree[u].push_back(newNode);
                    tree[newNode].push_back(u);

                    q.push(make_pair(newNode,make_pair(prev.first,i-1)));
                }
                else // add 2 neighbors
                {
                    tree.push_back(empVec);
                    int newNode=tree.size()-1;
                    tree[u].push_back(newNode);
                    tree[newNode].push_back(u);

                    q.push(make_pair(newNode,make_pair(prev.first,i-1)));


                    tree.push_back(empVec);
                    newNode=tree.size()-1;
                    tree[u].push_back(newNode);
                    tree[newNode].push_back(u);

                    q.push(make_pair(newNode,make_pair(prev.second,i-1)));
                }
            }
        }

        outputGraph(tree);
    }
    
    cout << "Finished output for case where the tree has a single centre." << endl;

    /*for(int i=1; i<=20; i++)
    {
        for(int j=1; j<=10; j++)
        {
            cerr << "dW[" << i << "][" << j << "]: " << dW[i][j] << endl;
        }
    }*/
    //cerr << dW[208][16] << " " << dW[209][16] << " " << dW[210][16] << endl;
    return 0;
}
