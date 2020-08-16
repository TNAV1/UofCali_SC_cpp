#include <iostream> // for cout cin
#include <ctime> // for rand seeding
#include <cstdlib> // rand
#include <vector> // for neighbour reporting
#include <iomanip> // for setprecision & fixed
#include <iterator> // for file loading
#include <fstream>  // for file loading
//#include <array> // for determining size

using namespace std; // for cout, endl

class Graph{
	
  public:
  int size; 
  double** matrix; // pointer to a an array of pointers (arrarys are automatically pointers)
  int txtsize;

  // Constructor to bring in the size variable
  Graph(int x) {size = x;}

  // create the matrix
  void create(double d, double mc, bool undirected)
  {
    matrix = new double*[size]; // pointer to heap memory of a certain number of arrays
    srand(time(0));
    double density = d;
    int q =0;
    for(int p = 0; p < size; ++p) {
      matrix[p] = new double[size]; // heap memory for an array of a certain size
    }
    for(int i = 0; i < size; ++i) { // only populate half the matrix if undirected
      if (undirected) {
        q=i;
      }else{
        q=0	;
      }
      for(int j = q; j < size; ++j) {
        // diagonal = false
        if(i == j) {
          matrix[i][j] = 0.0;
        }
        // symmetric
        else {
          double prob = rand() % 100/100.0;
          double cost = (rand() %100/100.0*(mc-1))+1; // rand cost 1.0 to maxcost (mc)
          if(prob < density){
            if (undirected) {
              matrix[i][j] = matrix[j][i] = cost;
            }else{
              matrix[i][j]  = cost;
            }  
          }else{
            if (undirected){
              matrix[i][j] = matrix[j][i] = 0.0;
            }else{
              matrix[i][j]  =  0.0;
            }  
          }    
        }
      }
    }
  }
  // calculate the density of the matrix
  double density(){
  	int edgecnt = 0;
    for(int i = 0; i < size; ++i) {
      for(int j = 0; j < size; ++j) {
      	if (matrix[i][j] !=0){
      		edgecnt=edgecnt+1;
      	}
      }
    } 
    double maxedge = (size-1)*size;
    double d = edgecnt/maxedge;
    return d;
  }
  
  // display the matrix
  void display(){
    for(int i = 0; i < size; ++i) {
      for(int j = 0; j < size; ++j) {
      	cout<<setprecision(1)<<fixed; // show 0.0 in the graph output instead of 0
        cout << matrix[i][j] << "  ";
      }
      cout << endl;
    }
  }
  
  // number of vertices in the graph
  int V() {return size;}

  // number of edges in the graph
  int E() {
  	  int edges = 0; 
  	  for(int i = 0; i < size; ++i) {
        for(int j = 0; j < size; ++j) {
        	if (matrix[i][j] != 0) {
        		edges = edges + 1;
        	}
        }
  	  }  
  	return edges;
  }
  
  // tests whether there is an edge from node x to node y
  bool adjacent(int i, int j) {
  	int adj = 0;
  	if (matrix[i][j] != 0){
  		adj = 1 ; 
  	}
  	return adj;
  }
  
  // lists all of the verticies connected to a given vertice
  vector <int> neighbour(int v){
  	vector <int> vctr;
  	for (int j=0; j<size; ++j){
    	if (matrix[v][j] !=0){
    	  vctr.push_back(j);  		
    	}
  	}
  	return vctr;
  	vctr.clear();
  }
  
  // add a vertice manually
  void add(int i, int j, double k){
  	if (matrix[i][j] != 0){
  	 cout << "WARNING -- [" << i << "," << j << "] connection already exists" << endl;
  	} else{
     matrix[i][j] = k;
  	}
  }
  
  // delete a vertice manually
  void VertDel (int m, int n){
  	if (matrix[m][n] == 0){
  	 cout << "WARNING -- [" << m << "," << n << "] never existed" << endl;
  	} else{
     matrix[m][n] = 0.0;
  	}
  }
  
  // get edge value
  double get_edge_value(int i, int j){return matrix[i][j];}
  
  // set edge value (same as add but will ovewrite a current cost value)
  void set_edge_value(int i, int j, double k){matrix[i][j] = k;}
  
  
  double ** getGraph(){
         return matrix;
     }
  
  double AveGraphLen(){
  	double sumcost=0;
  	double cntcost=0;
  	double cost=0;
    for (int i=0; i<size; ++i)	{
      for (int j=0; j<size; ++j)	{
      	cost = matrix[i][j];
      	if (cost!=0){
      		sumcost=sumcost+cost;
      		cntcost=cntcost+1;
      	}
      }
    }
    cout << "AverGraphLen " << sumcost/cntcost << endl;
    return sumcost/cntcost;
  }

  double ** LoadMatrix(){

    ifstream data_file("w:\\CPP\\LoadFile\\data.txt");
    istream_iterator<int> start(data_file), end;
    vector<int> data(start,end);
	
    int ds = data.size();
    if (data.size() == 0){
      cout << "ERROR" << endl;
      cout << "file does not exist" << endl;
    } else {

      // create an empty matrix to hold all the data
      txtsize = data[0];
      matrix = new double*[txtsize];
      for(int p = 0; p < txtsize; ++p) {
        matrix[p] = new double[txtsize]; // heap memory for an array of a certain txtsize
      }
      for(int p = 0; p < txtsize; ++p) { // init the matrix to zeros
        for(int q = 0; q < txtsize; ++q) {
          matrix[p][q] = 0; 
        }
      }
      for(int i=1; i<=ds-3; i=i+3){
        matrix[data[i]][data[i+1]]=data[i+2];
      }
    }
    return matrix;
    }
  
  // free the heap memory
  void del(){
	  delete []matrix;
	  cout << "Graph heap memory cleared" << endl;
  }
  
};

// override << to allow printing of a vector
template <typename T> 
ostream& operator<<(ostream& os, const vector<T>& vec) 
{ 
    os << "["; 
    for (int i = 0; i < vec.size(); ++i) { 
        os << vec[i]; 
        if (i != vec.size() - 1) 
            os << ", "; 
    } 
    os << "]\n"; 
    return os; 
} 



class ShortestPath{
  public:
  
    // Constructor to bring in the matrix/graph to be evaluated
    double ** mtrx; // use pointer to large matrix to avoid using extra memory
    vector <int> OpenBin;
    vector <int> ClosedBin;
    vector <double> MinDist;
    vector <int> Parent;
	int msize; 
	int bwsize; // backwards path size
    
	ShortestPath(double ** x, int y){ // use pointer to large matrix to avoid using extra memory
		mtrx = x;
		msize = y;
	}
	
	int NumCol(){
//		int len = mtrx.length;
 //     int len=size(mtrx[0]);
        int len = sizeof( *mtrx ) / sizeof( **mtrx );
//     	cout << "length " << len << endl;
		return len;
	}
	
	
	// lists all of the verticies connected to a given vertice
    vector <int> neighbour(int v){
  	  vector <int> spvctr;
  	  NumCol();
  	  for (int j=0; j<msize; ++j){ //TODO: measure the size of the matrix and use the column size
      	if (mtrx[v][j] !=0){
      	  spvctr.push_back(j);  		
      	}
  	  }
  	  return spvctr;
  	  spvctr.clear();
    }
    
    
    //Check if vertice is in the bin yet
    bool InBin(int v, vector <int> Bin) {
        bool found = 0;
        for (int j=0; j < Bin.size(); ++j) {
          if (v==Bin[j]) {
            found = 1;
            break; 
          }
        }
        return found;
        Bin.clear();
     }
  
  
    //Calculate the shortest path
    void CalcPath(int u, int w){
      vector <int> conn;
      int cnt = 0;
      int loc =0;
      int parentloc=0;
      double md; 
      int mdvert;
      double CalcDist;
      bool found = 0;
      
      int CurrVert = u; // set the start value 
      OpenBin.push_back(u); // add the start vertex to the open bin
      ClosedBin.push_back(u); // the start vertex is not closed yet
      MinDist.push_back(0); // the start vertex gets a distance of zero
      Parent.push_back(999999); // the start vertex has no parent
               
      while (!InBin(w,ClosedBin) && cnt<50){ // remove safety cnt 
        conn = neighbour(CurrVert); // bring in the connected vertices
        for (int i=0; i < conn.size(); ++i ){ // cycle thru all connected vertices
          if (!InBin(conn[i],ClosedBin)){ // do the following steps only if it is not in the ClosedBin
            found = InBin(conn[i],OpenBin); 
            if (found!=1){ //if not found in the OpenBin or ClosedBin make a spot for it in all 4 tracking vecotrs
            	OpenBin.push_back(conn[i]); 
            	ClosedBin.push_back(999999); 
            	MinDist.push_back(999999); 
            	Parent.push_back(999999); 
            }
            loc = VectorFind(OpenBin, conn[i]); // find the neighbour location in the OpenBin
            parentloc = VectorFind(OpenBin, CurrVert); // find the parent location in the OpenBin
            CalcDist = mtrx[CurrVert][conn[i]] + MinDist[parentloc];
            if (CalcDist < MinDist[loc]){
            	MinDist[loc] = CalcDist; //set the min distance
            	Parent[loc] = CurrVert;  // set the parent
            }
          }
        }
        parentloc = VectorFind(OpenBin, CurrVert); // find the parent location in the OpenBin
        ClosedBin[parentloc] = CurrVert; // set vertex as closed 
       // cout << CurrVert << " Has these neighbours " << conn << endl;
       // cout << "OpenBin   " << OpenBin << endl;
       // cout << "ClosedBin " << ClosedBin << endl;
       // cout << "MinDist   " << MinDist << endl;
       // cout << "Parent    " << Parent << endl;
        
        // select the smallest distance that is not closed to evalaute next
        md = 999999;
        mdvert=0;
        for (int j=1; j < MinDist.size(); ++j){ // avoid the start vertex which as distance zero
          	if (MinDist[j]<md) {
          		if (ClosedBin[j]!=OpenBin[j]) { // dont use it if it is in the Closed bin
          		  md=MinDist[j];
          		  mdvert=OpenBin[j];
          		}
          	}
        }
       // cout << "Next min dist " << md << " " << mdvert << endl;
        CurrVert = mdvert;
        cnt=cnt+1; // remove safety cnt
      }
      conn.clear();
  }
  
  
  double path_size(int u, int w){
    double MinCost;
    // find the end vertex location in the OpenBin
    int loc = VectorFind(OpenBin, w);
    MinCost = MinDist[loc];
  	return MinCost;
  }
  
  bool connected(int u, int w){
  	  double MinCost;
  	  bool connect;
      MinCost = path_size(u,w);
      if (MinCost == 0){
        connect = 0; 
      }else{
        connect = 1;
      }
  	return connect;
  }
  
  
  vector <int> path(int u, int w){
  	  vector <int> spathbw; // shortest path listed from end to start
  	  int NextBWVert; // next vertice working backwards
  	 // int bwsize; // backwards path size
  	  vector <int> spath;  // shortest path listed from start to end
  	  int vert; 
  	  
      spathbw.push_back(w); // build the path working backwards strating from the end vertex
      // find the parent
      NextBWVert=w;
      int pcntr = 0;
      while (NextBWVert != u && pcntr<msize ){
        int loc = VectorFind(OpenBin, NextBWVert); // find the end vertex location in the OpenBin
        NextBWVert=Parent[loc]; // look up the parent at that same location
        spathbw.push_back(NextBWVert); // add it to the backwards path list
        pcntr=pcntr+1; // counter to prevent getting stuck in the while loop
      }
      // create the path from start to end
      bwsize=spathbw.size();
      for (int j=0; j < bwsize; ++j){ spath.push_back(spathbw[bwsize-1-j]);}
      //cout << "shortest path is " << spath << endl;
      return spath;
      spathbw.clear();
      spath.clear();
  }
  

  int VectorFind(vector <int> vec, int x){
    // find the x location in the vector
    int loc = 0;
    for (int j=0; j < vec.size(); ++j){
      if (vec[j]==x) {loc = j;}
    }
    return loc;
    vec.clear();
  }
  
  // free the heap memory
  void del(){
	  //delete []mtrx;
	  OpenBin.clear();
      ClosedBin.clear();
      MinDist.clear();
      Parent.clear();
	  cout << "ShortestPath heap memory cleared" << endl;
  }

};


class MinSpanTree{
  public:
    // Constructor to bring in the matrix/graph to be evaluated
    double ** mtrx; // use pointer to large matrix to avoid using extra memory
    int msize; // size of matrix
    vector <int> OpenBin;
    vector <int> ClosedBin;
    vector <double> MinDist;
    vector <int> Parent;
    
	MinSpanTree(double ** x, int y){ // use pointer to large matrix to avoid using extra memory
		mtrx = x;
		msize = y;
	}
	
	// lists all of the verticies connected to a given vertice
    vector <int> neighbour(int v){
  	  vector <int> spvctr;
  	  for (int j=0; j<msize; ++j){ //TODO: measure the size of the matrix and use the column size
      	if (mtrx[v][j] !=0){
      	  spvctr.push_back(j);  		
      	}
  	  }
  	  return spvctr;
  	  spvctr.clear();
    }
    
    //Check if vertice is in the bin yet
    bool InBin(int v, vector <int> Bin) {
        bool found = 0;
        for (int j=0; j < Bin.size(); ++j) {
          if (v==Bin[j]) {
            found = 1;
            break; 
          }
        }
        return found;
        Bin.clear();
     }
	
	void MSTcalc(double ** m){
	  vector <int> conn;
	  double MinConn; 
      int CurrVert=0;
      bool done = 0;
      int SmallDist;
      int loc;
	  
      OpenBin.push_back(0); // add the start vertex to the open bin
      ClosedBin.push_back(0); // the start vertex is not closed yet
      MinDist.push_back(0); // the start vertex gets a distance of zero
      Parent.push_back(999999); // the start vertex has no parent
      
      while (!done){
  	  	conn = neighbour(CurrVert); //find the neighbours of the vertex
  	  	for (int j=0; j < conn.size(); ++j){
  	  	  if (m[CurrVert][conn[j]]!=0){ // Only process lenghts that are not zero
            if(!InBin(conn[j],OpenBin) && !InBin(conn[j],ClosedBin)){//If not in a bin added the neighbour
  	  	      OpenBin.push_back(conn[j]);
  	  	      Parent.push_back(CurrVert);
  	  	      MinDist.push_back(m[CurrVert][conn[j]]);
  	  	      ClosedBin.push_back(999999);
  	  	    } else if (InBin(conn[j],OpenBin) && !InBin(conn[j],ClosedBin) && MinDist[VectorFind(OpenBin,conn[j])]>m[CurrVert][conn[j]] ){ //If already in OpenBin but new edge has a smaller value
  	  	      Parent[VectorFind(OpenBin,conn[j])]=CurrVert;
  	  	      MinDist[VectorFind(OpenBin,conn[j])]=m[CurrVert][conn[j]];
            }
  	  	  }
  	  	  ClosedBin[VectorFind(OpenBin,CurrVert)]=CurrVert; // close the parent vertex
  	  	}  	  	
  	  	// pick the small distance that is not closed
	    SmallDist = 999999;
	    for (int k=0; k < MinDist.size(); ++k){
          if (MinDist[k]<SmallDist && OpenBin[k]!=ClosedBin[k]){ //if it is a min distance and not closed
          	SmallDist = MinDist[k];
          	loc = k;
          }
	    }
	    if (OpenBin==ClosedBin){
	    	done=1;
	    } // if there is nothing more to evalute stop the loop
	    CurrVert=OpenBin[loc];
      }    

      conn.clear();
	} 
	
	double MSTvalue(){
	  double TreeValue=0.0;
      for (int i=1; i<MinDist.size(); ++i){
      	TreeValue=TreeValue+MinDist[i];
      }
      return TreeValue;
	}
	
	void MSTEdges(){
	  cout << "Min Spanning Tree edge list is ... " ;
	  for (int i=1; i<Parent.size(); ++i){
	  	cout << "(" << Parent[i] << "," << ClosedBin[i] << ") ";
	  }
	  cout << endl;
	}
	
    int VectorFind(vector <int> vec, int x){
      // find the x location in the vector
      int loc = 0;
      for (int j=0; j < vec.size(); ++j){
        if (vec[j]==x) {loc = j;}
      }
      return loc;
      vec.clear();
    }
  
  // free the heap memory
  void del(){
	  OpenBin.clear(); 
    ClosedBin.clear();
    MinDist.clear();
    Parent.clear();
	  cout << "MST heap memory cleared" << endl;
  }
};
  
// *******************************************************
//
// Main
//
// *******************************************************

int main()
{
  int size = 4;
  //cout << "Enter the number of vertices " << endl;
  //cin >> size;
  //cout << "You entered " << size << endl;
  Graph G(size);       // declaration of type Graph and variable G
  
  // create the graph with appropriate edge density & cost
  double dens = 0.8;
  double maxcost = 10;
  bool undirected = 1; // undirected means same cost in both directions
  G.create(dens, maxcost, undirected);
  cout<<setprecision(1)<<fixed; // show 0.0 in the graph output instead of 0
  cout << "Calculated density is " << G.density() << " vs requested density of " << dens << endl;
  G.display();
  
  G.AveGraphLen();
  
  
  // test for number of vertices & edges
  //cout << "Number of vertices " << G.V() << endl;
  //cout << "Number of edges " << G.E() << endl;
  
  // test for connected verticies
  int i = 0; 
  int j = 0;
  //cout << "Enter two veritices to see if they are connected" << endl;
  //cout << "Vertice 1" << endl;
  //cin >> i;
  //cout << "Vertice 2" << endl;
  //cin >> j;
  //cout << "Vertices " << i << " and " << j << " are connected " << G.adjacent(i,j)  << endl;
  
  // test for neighbouring verticies
  int testvertice = 1;
  //cout << "Enter a vertice to see a list of neighbours that you could travel to" << endl;
  //cin >> testvertice;
  //cout << "Vertice " << testvertice << " is connected to vertices " << G.neighbour(testvertice) << endl;

  // add a vertice manually
  int va1 = 0;
  int va2 = 0;
  double vvalue = 4.0;
  //cout << "Enter 2 vertices to add a connection between" << endl;
  //cin >> va1;
  //cin >> va2;
  //cout << "Enter the cost value for that edge" << endl;
  //cin >> vvalue;
  //cout << "Cost =" << vvalue << endl;
  if (va1<size && va2<size){ // prevent user from enter an unrealistic vertices
    G.add(va1,va2,vvalue);
    //cout << "Current matrix is..." << endl;
    //G.display();   	
  } else {
  	//cout << "ERROR - vertice " << va1 << " or " << va2 << " is larger than the max possible value of " << size-1 << endl;
  }

  
  // delete a vertice manually
  int vd1 = 0;
  int vd2 = 0;
  //cout << "Enter 2 vertices to delete a connection between" << endl;
  //cin >> vd1;
  //cin >> vd2;
  if (vd1<size && vd2<size){ // prevent user from enter an unrealistic vertices
    G.VertDel(vd1,vd2);
    //cout << "Current matrix is..." << endl;
    //G.display();
  } else {
    //cout << "ERROR - vertice " << vd1 << " or " << vd2 << " is larger than the max possible value of " << size-1 << endl;
  }
  
  // get an edge value
  int ev1 = 0;
  int ev2 = 0;
  //cout << "Select an edge by defining 2 verticies, its value will be reutrned " << endl;
  //cin >> ev1;
  //cin >> ev2;
  if (ev1<size && ev2<size){ // prevent user from enter an unrealistic vertices
    //cout << "Edge value is " << G.get_edge_value(ev1, ev2) << endl;
  } else {
    //cout << "ERROR - vertice " << ev1 << " or " << ev2 << " is larger than the max possible value of " << size-1 << endl;
  }
  
  // set an edge value
  int evs1 = 0;
  int evs2 = 0;
  double setev = 0; 
  //cout << "Select an edge to update its cost value by defining 2 verticies" << endl;
  //cin >> evs1;
  //cin >> evs2;
  //cout << "Enter the cost value " << endl;
  //cin >> setev;
  if (evs1<size && evs2<size){ // prevent user from enter an unrealistic vertices
    G.set_edge_value(evs1,evs2,setev);
    //cout << "Current matrix is..." << endl;
    //G.display();
  } else {
    //cout << "ERROR - vertice " << evs1 << " or " << evs2 << " is larger than the max possible value of " << size-1 << endl;
  }
  
  // bring the graph into main as a pointer so it does not consume excessive memory
  double ** graph =   G.getGraph();
  //cout << "Main graph: " << graph[1][3] << endl;
  
  // pass the matrix to be evaluated into the class ShortestPath
  ShortestPath sp(graph, size);
  
  // Determine the shortest path
 
  int vsp1 = 0;
  int vsp2 = size-1;
  cout << "Finding shortest path from " << vsp1 << " to " << vsp2 << "....." << endl;
  sp.CalcPath(vsp1,vsp2);
  if (!sp.connected(vsp1,vsp2)){
    cout << "ERROR - " << vsp1 << " and " << vsp2 << " are not connected !!!!" << endl;
  }else{
    cout << "Min cost of path is " << sp.path_size(vsp1,vsp2) << endl;
    cout << "shortest path is " << sp.path(vsp1,vsp2) << endl;
    cout << "Path lenght is " << sp.bwsize << endl;
  }
  
   double ** txtmtrx = G.LoadMatrix();
   int txtsize = G.txtsize;

    for(int p = 0; p < 20; ++p) {    // print the data file matrix
      for(int q = 0; q < 20; ++q) {
        cout << txtmtrx[p][q] << " " ;
      }
      cout << endl;
    }  

   MinSpanTree MST(txtmtrx, txtsize);
   MST.MSTcalc(txtmtrx);
   cout << "Minimum spanning tree lenght is " << MST.MSTvalue() << endl;
   MST.MSTEdges();

  // clear the heap
  G.del(); 
  sp.del();
  MST.del();
  
  return 0;
}



