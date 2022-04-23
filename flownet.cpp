#include "Rcpp.h"
#include <vector>
#include <string>
#include <map>

// TODO : SOLO PARA CICLOS : NO REVISAR DESDE UN NODO INICIAL NODOS MAS CHICOS (XQ TAMBIEN SE REVISARAN POR EL OTRO LADO)
// TODO : MANTENER UN MAPA DE SETS DE CICLOS POR NODO INICIAL, QUIZAS HASTA SEGUNDO ORDEN INCLUSO, PARA DISMINUIR LAS INSERCIONES
// TODO : HACER UNA FUNCION APARTE DE TRAYECTORIAS?

using namespace std;
using namespace Rcpp; 

namespace std {
template <>
struct hash<vector<int> >
{
  typedef vector<int>      argument_type;
  typedef std::size_t  result_type;
  
  result_type operator()(const vector<int> & t) const
  {
    std::size_t seed = 0;//t.size();
    for(auto& i : t) {
      //seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      seed += pow(2,i);
    }
    return seed;
  }
};
}

map<int,vector<int> > to_elements(NumericMatrix m){
  map<int,vector<int> > out;
  for(int i=0;i<m.rows();++i){
    vector<int> current;
    for(int j=0;j<m.cols();++j){
      if(j==i)
        continue;
      if(m(i,j)>0){
        current.push_back(j);
      }
    }
    out[i]=current;
  }
  return out;
}
void add_cycle(unordered_set<vector<int> >* cumulative_cycles,vector<int> new_cycle){
  int min = 100000;
  int index=0;
  for(int i=0;i<new_cycle.size();++i){
    if(min>new_cycle[i]){
      min=new_cycle[i];
      index=i;
    }
  }
  vector<int> cycle;
  int c_size = new_cycle.size();
  for(int i=0;i<c_size;++i){
    cycle.push_back(new_cycle[(index+i)%c_size]);
  }
  int pre = cumulative_cycles->size();
  cumulative_cycles->insert(cycle);
  if(pre<cumulative_cycles->size() && (cumulative_cycles->size()%10000)==0){
    cout << "CYCLES AT : " << cumulative_cycles->size() << endl;
  }
}
vector<vector<int> > next_paths(map<int,vector<int> > connected, 
                                vector<vector<int> > current_paths,
                                vector<vector<int> >* cumulative_trajectories,
                                unordered_set<vector<int> >* cumulative_cycles,
                                int output_node){
  vector<vector<int> > out;
  for (auto p = begin (current_paths); p != end (current_paths); ++p) {
    vector<int> which_connected = connected[p->back()];
    for(auto c = begin(which_connected);c!=end(which_connected);++c){
      if((*c)==output_node){
        vector<int> done = *p;
        done.push_back(*c);
        for(auto t=begin(done);t!=end(done);++t){
          (*t) = (*t) + 1;
        }
        cumulative_trajectories->push_back(done);
      }else{
        bool hasit=false;
        for (auto n = begin (*p); n != end (*p); ++n){
          if((*n)==(*c)){
            hasit=true;
            vector<int> cycle;
            for (auto cy = n; cy != end (*p); ++cy){
              cycle.push_back((*cy)+1);
            }
            add_cycle(cumulative_cycles,cycle);
            break;
          }
        }
        if(!hasit){
          vector<int> done = *p;
          done.push_back(*c);
          out.push_back(done);
        }
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
vector<vector< vector<int> > > paths(NumericMatrix m,int input_node,int output_node){
  vector<vector<int> > out;
  unordered_set< vector<int> > cycles; //
  //cout << m << endl;
  map<int,vector<int> > connected = to_elements(m);
  vector<vector<int> > paths;
  vector<int> conn = connected[input_node-1];
  for(auto it=begin(conn);it!=end(conn);++it){
    vector<int> p;
    p.push_back(input_node-1);
    p.push_back(*it);
    if(*it == (output_node-1)){
      out.push_back(p);
    }else{
      paths.push_back(p);
    }
  }
  while(paths.size()>0){
    paths = next_paths(connected,paths,&out,&cycles,output_node-1);
  }
  vector<vector< vector<int> > > final_out;
  final_out.push_back(out);
  vector< vector<int> > cycles_out;
  for(auto cycle : cycles){
    vector<int> c;
    for(auto n=begin(cycle);n!=end(cycle);n++){
      c.push_back(*n);
    }
    cycles_out.push_back(c);
  }
  final_out.push_back(cycles_out);
  return final_out;
}

/*** R
system.time(res<-paths(matrix(sample(c(0,1),100,replace = TRUE),nrow=10,ncol=10),9,10))
print(length(res))
*/
