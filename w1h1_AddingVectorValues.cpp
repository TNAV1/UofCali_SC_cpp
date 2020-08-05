// Convert this program to C++
// change to C++ io
// change to one line comments
// change defines of constants to const
// change array to vector<>
// inline any short function

#include <iostream> // supports cout commands
#include <vector> // supports use of vectors
using namespace std; // shortens command names

const int N = 40 ; // convert from #define to const

template <class T> // template used so array or vector could be passed
inline void sum(int*p, int n, T d) { // T defines the type used for d	
  *p = 0;
  for(int i = 0; i < n; ++i) // local declaration of i for tighter code
    *p = *p + d[i];
}

int main()
{
  int accum = 0;
  vector <int> data; //convert array to vector
  for(int i = 0; i < N; ++i) // local declaration of i for tighter code
   data.push_back(i); //build vector
   sum(&accum, N, data);
   cout << "sum is " << accum << endl; //use cout instead of printf
  return 0;
}
