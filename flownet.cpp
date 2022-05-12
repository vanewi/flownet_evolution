#include "Rcpp.h"
#include <vector>
#include <string>
#include <map>

// TODO : SOLO PARA CICLOS : NO REVISAR DESDE UN NODO INICIAL NODOS MAS CHICOS (XQ TAMBIEN SE REVISARAN POR EL OTRO LADO)
// TODO : MANTENER UN MAPA DE SETS DE CICLOS POR NODO INICIAL, QUIZAS HASTA SEGUNDO ORDEN INCLUSO, PARA DISMINUIR LAS INSERCIONES
// TODO : HACER UNA FUNCION APARTE DE TRAYECTORIAS?
// TODO : COMPRESS CYCLES USING TWO HASHES: ONE FOR WHICH NODES AND ONE FOR ORDER
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
    std::size_t seed = 0;
    for(auto& i : t) {
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
void add_cycle(map<int,unordered_set< vector<int> > >* cumulative_cycles,vector<int> new_cycle,NumericMatrix* agg_cycle){
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
  if(cumulative_cycles->find(min)==cumulative_cycles->end()){
    (*cumulative_cycles)[min]=unordered_set< vector<int> >();
  }
  int pre = (*cumulative_cycles)[min].size();
  (*cumulative_cycles)[min].insert(cycle);
  if(pre<(*cumulative_cycles)[min].size() ){
    if(((*cumulative_cycles)[min].size()%10000)==0){
      int sum=0;
      for(auto mit=cumulative_cycles->begin();mit!=cumulative_cycles->end();mit++){
        sum+=mit->second.size();
      }
      cout << "CYCLES AT : " << sum << endl;
    }
    for(int i=1;i<c_size;++i){
      (*agg_cycle)(cycle[i-1]-1,cycle[i]-1)=(*agg_cycle)(cycle[i-1]-1,cycle[i]-1)+1;
    }
    (*agg_cycle)(cycle.back()-1,cycle[0]-1)=(*agg_cycle)(cycle.back()-1,cycle[0]-1)+1;
  }
}
vector<vector<int> > next_paths(map<int,vector<int> > connected, 
                                vector<vector<int> > current_paths,
                                vector<vector<int> >* cumulative_trajectories,
                                map<int,unordered_set<vector<int> > >* cumulative_cycles,
                                int output_node,
                                NumericMatrix* agg_cycle,
                                NumericMatrix* agg_traj,
                                bool save_trajectories){
  vector<vector<int> > out;
  for (auto p = begin (current_paths); p != end (current_paths); ++p) {
    vector<int> which_connected = connected[p->back()];
    for(auto c = begin(which_connected);c!=end(which_connected);++c){
      if((*c)==output_node){
        vector<int> done = *p;
        done.push_back(*c);  
        int pre_node;
        bool started=false;
        for(auto t=begin(done);t!=end(done);++t){
          if(!started){
            started=true;
            pre_node=*t;
          }
          else{
            (*agg_traj)(pre_node,*t)=(*agg_traj)(pre_node,*t)+1;
            pre_node = *t;
          }
          (*t) = (*t) + 1;
        }
        if(save_trajectories){
          cumulative_trajectories->push_back(done);
        }else{
          if(cumulative_trajectories->size()==0){
            vector<int> acum;
            acum.push_back(1);
            cumulative_trajectories->push_back(acum);
          }else{
            cumulative_trajectories->back()[0]=cumulative_trajectories->back()[0]+1;
          }
        }
      }else{
        bool hasit=false;
        for (auto n = begin (*p); n != end (*p); ++n){
          if((*n)==(*c)){
            hasit=true;
            vector<int> cycle;
            for (auto cy = n; cy != end (*p); ++cy){
              cycle.push_back((*cy)+1);
            }
            add_cycle(cumulative_cycles,cycle,agg_cycle);
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
// Buscar trayectorias y ciclos desde el nuevo link
// Si en la busqueda encontramos un ciclo que no contenga el link, se desecha
// El resultado de los ciclos nuevos se une a los que teniamos
// Las trayectorias nuevas solo aportan "multiplicando' las antiuguas 'unicas' truncadas al nodo inicial, que intersectan solamente a ese nodo inicial del link nuevo
// Si sacamos un link se eliminan todas las trayectorias que continenen ese link y los ciclos que contienen ese link



//METRICAS CICLOS
// Participacion promedio por link de ciclos (media de numero de ciclos que participa un link sobre el total)
// Media de participacion por nodo de ciclos
// Identificacion de indegrees y outdegrees por ciclos (solo desde y hacia nodos fuera del ciclo)

// [[Rcpp::export]]
vector<vector< vector<int> > > paths(NumericMatrix m,int input_node,int output_node,bool save_trajectories=true){
  NumericMatrix aggregated_cycle_participation( m.rows() );
  NumericMatrix aggregated_trajectory_participation( m.rows() );
  vector<vector<int> > trajectories;
  map<int,unordered_set< vector<int> > > cycles; //
  //cout << m << endl;
  map<int,vector<int> > connected = to_elements(m);
  vector<vector<int> > paths;
  vector<int> conn = connected[input_node-1];
  for(auto it=begin(conn);it!=end(conn);++it){
    vector<int> p;
    p.push_back(input_node-1);
    p.push_back(*it);
    if(*it == (output_node-1)){
      trajectories.push_back(p);
    }else{
      paths.push_back(p);
    }
  }
  while(paths.size()>0){
    paths = next_paths(connected,paths,&trajectories,&cycles,output_node-1,&aggregated_cycle_participation,&aggregated_trajectory_participation,save_trajectories);
  }
  vector<vector< vector<int> > > final_out;
  vector< vector<int> > cycles_out;
  for(auto mit=cycles.begin();mit!=cycles.end();++mit){
    for(auto cycle : mit->second){
      vector<int> c;
      for(auto n=begin(cycle);n!=end(cycle);n++){
        c.push_back(*n);
      }
      cycles_out.push_back(c);
    }
  }
  vector<vector<int> > agg_cycles,agg_trajs;
  for(int i=0;i<aggregated_cycle_participation.rows();++i){
    vector<int> cycle_row,traj_row;
    for(int j=0;j<aggregated_cycle_participation.cols();++j){
      cycle_row.push_back(aggregated_cycle_participation(i,j));
      traj_row.push_back(aggregated_trajectory_participation(i,j));
    }
    agg_cycles.push_back(cycle_row);
    agg_trajs.push_back(traj_row);
  }
  final_out.push_back(trajectories); //res[[1]]
  final_out.push_back(agg_trajs); //res[[2]]
  final_out.push_back(cycles_out); //res[[3]]
  final_out.push_back(agg_cycles); //res[[4]]
  return final_out;
}

/*** R
system.time(res<-paths(matrix(sample(c(0,1),100,replace = TRUE),nrow=10,ncol=10),9,10))
print(length(res))
*/
