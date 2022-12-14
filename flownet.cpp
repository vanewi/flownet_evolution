#include "Rcpp.h"
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <vector>
#include <string>
#include <map>

// TODO : SOLO PARA CICLOS : NO REVISAR DESDE UN NODO INICIAL NODOS MAS CHICOS (XQ TAMBIEN SE REVISARAN POR EL OTRO LADO)
// TODO : MANTENER UN MAPA DE SETS DE CICLOS POR NODO INICIAL, QUIZAS HASTA SEGUNDO ORDEN INCLUSO, PARA DISMINUIR LAS INSERCIONES
// TODO : HACER UNA FUNCION APARTE DE TRAYECTORIAS?
// TODO : COMPRESS CYCLES USING TWO HASHES: ONE FOR WHICH NODES AND ONE FOR ORDER
// TODO : MAKE HASH DURING PATH IN ORDER NOT TO ITERATE THROUGH EVERY ELEMENT IN CURRENT PATH TO CHECK CYCLE
using namespace std;
using namespace Rcpp; 
using namespace RcppParallel;

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
//    Length   first    second       cycle
//map<int, map<int, map<int, unordered_set< vector<int> > > > >
void add_cycle(map<int, map<int, map<int, unordered_set< vector<int> > > > >* cumulative_cycles,vector<int> new_cycle,NumericMatrix* agg_cycle,int* totCycles){
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
  int sec = cycle[1];
  int pre = (*cumulative_cycles)[c_size][min][sec].size();
  (*cumulative_cycles)[c_size][min][sec].insert(cycle);
  if(pre<(*cumulative_cycles)[c_size][min][sec].size() ){
    *totCycles = *totCycles + 1;
    if(*totCycles%10000==0){
      cout << "CYCLES AT : " << *totCycles << endl;
    }
    for(int i=1;i<c_size;++i){
      (*agg_cycle)(cycle[i-1]-1,cycle[i]-1)=(*agg_cycle)(cycle[i-1]-1,cycle[i]-1)+1;
    }
    (*agg_cycle)(cycle.back()-1,cycle[0]-1)=(*agg_cycle)(cycle.back()-1,cycle[0]-1)+1;
  }
}

void add_cycle_parallel(tthread::mutex* cycle_mutex,map<int, map<int, map<int, unordered_set< vector<int> > > > >* cumulative_cycles,vector<int> new_cycle,NumericMatrix* agg_cycle,int* totCycles){
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
  int sec = cycle[1];
  cycle_mutex->lock();
  int pre = (*cumulative_cycles)[c_size][min][sec].size();
  (*cumulative_cycles)[c_size][min][sec].insert(cycle);
  if(pre<(*cumulative_cycles)[c_size][min][sec].size() ){
    *totCycles = *totCycles + 1;
    if(*totCycles%10000==0){
      cout << "CYCLES AT : " << *totCycles << endl;
    }
    for(int i=1;i<c_size;++i){
      (*agg_cycle)(cycle[i-1]-1,cycle[i]-1)=(*agg_cycle)(cycle[i-1]-1,cycle[i]-1)+1;
    }
    (*agg_cycle)(cycle.back()-1,cycle[0]-1)=(*agg_cycle)(cycle.back()-1,cycle[0]-1)+1;
  }
  cycle_mutex->unlock();
}

struct NextPathsWorker : public Worker
{
  // source matrix
  tthread::mutex* cycle_mutex;
  tthread::mutex* trajectory_mutex;
  tthread::mutex* path_mutex;
  
  // destination matrix
  vector<vector<int> > outPaths;
  map<int,vector<int> >* connected;
  vector<vector<int> >* current_paths;
  vector<vector<int> >* cumulative_trajectories;
  map<int, map<int, map<int, unordered_set< vector<int> > > > >* cumulative_cycles;
  int output_node;
  NumericMatrix* agg_cycle;
  NumericMatrix* agg_traj;
  bool save_trajectories;
  int* totCycles;
  
  void Restart(vector<vector<int> >* current_paths_){
    current_paths=current_paths_;
    outPaths.clear();
  }
  
  // initialize with source and destination
  NextPathsWorker(tthread::mutex* cycle_mutex,
                  tthread::mutex* trajectory_mutex,
                  tthread::mutex* path_mutex,
                  map<int,vector<int> >* connected, 
                  vector<vector<int> >* current_paths,
                  vector<vector<int> >* cumulative_trajectories,
                  map<int, map<int, map<int, unordered_set< vector<int> > > > >* cumulative_cycles,
                  int output_node,
                  NumericMatrix* agg_cycle,
                  NumericMatrix* agg_traj,
                  bool save_trajectories,int* totCycles) 
    : cycle_mutex(cycle_mutex), trajectory_mutex(trajectory_mutex),path_mutex(path_mutex),connected(connected),current_paths(current_paths),
      cumulative_trajectories(cumulative_trajectories),cumulative_cycles(cumulative_cycles),output_node(output_node),
      agg_cycle(agg_cycle),agg_traj(agg_traj),save_trajectories(save_trajectories),totCycles(totCycles){}
  
  // take the square root of the range of elements requested
  void operator()(std::size_t begin, std::size_t end) {
    //trajectory_mutex->lock();
    //cout << "PATH INTERVAL : " << (end-begin) << endl;
    //trajectory_mutex->unlock();
    vector<vector<int> > currentOut;
    for (auto p = current_paths->begin()+begin; p != current_paths->begin()+end; ++p) {
    vector<int> which_connected = (*connected)[p->back()];
    for(auto c = which_connected.begin();c!=which_connected.end();++c){
      if((*c)==output_node){
        /*vector<int> done = *p;
        done.push_back(*c);  
        int pre_node;
        bool started=false;
        for(auto t=done.begin();t!=done.end();++t){
          if(!started){
            started=true;
            pre_node=*t;
          }
          else{
            trajectory_mutex->lock();
            (*agg_traj)(pre_node,*t)=(*agg_traj)(pre_node,*t)+1;
            pre_node = *t;
            trajectory_mutex->unlock();
          }
          (*t) = (*t) + 1;
        }
        trajectory_mutex->lock();
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
        trajectory_mutex->unlock();*/
      }else{
        bool hasit=false;
        for (auto n = p->begin(); n != p->end(); ++n){
          if((*n)==(*c)){
            hasit=true;
            vector<int> cycle;
            for (auto cy = n; cy != p->end(); ++cy){
              cycle.push_back((*cy)+1);
            }
            add_cycle_parallel(cycle_mutex,cumulative_cycles,cycle,agg_cycle,totCycles);
            break;
          }
        }
        if(!hasit){
          vector<int> done = *p;
          done.push_back(*c);
          currentOut.push_back(done);
        }
      }
    }
  }
    path_mutex->lock();
    outPaths.insert(outPaths.end(), currentOut.begin(), currentOut.end());
    path_mutex->unlock();
  }
};

vector<vector<int> > next_paths(map<int,vector<int> >* connected, 
                                vector<vector<int> >* current_paths,
                                vector<vector<int> >* cumulative_trajectories,
                                map<int, map<int, map<int, unordered_set< vector<int> > > > >* cumulative_cycles,
                                int output_node,
                                NumericMatrix* agg_cycle,
                                NumericMatrix* agg_traj,
                                bool save_trajectories,int* totCycles){
  vector<vector<int> > out;
  for (auto p = current_paths->begin(); p != current_paths->end(); ++p) {
    vector<int> which_connected = (*connected)[p->back()];
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
            add_cycle(cumulative_cycles,cycle,agg_cycle,totCycles);
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
vector<vector< vector<int> > > paths(NumericMatrix m,int input_node,int output_node,bool save_trajectories=true,bool use_parallel=false,int limit_for_parallel=1000){
  NumericMatrix aggregated_cycle_participation( m.rows() );
  NumericMatrix aggregated_trajectory_participation( m.rows() );
  vector<vector<int> > trajectories;
  map<int, map<int, map<int, unordered_set< vector<int> > > > > cycles; //
  for(int i=2;i<(m.rows()-1);++i){// from 2 to nodes-2 (input and output)
    cycles[i]=map<int, map<int, unordered_set< vector<int> > > >();
    for(int j=1;j<(m.rows()-2);++j){// from 1 to core_nodes-1 (smallest to the one previous to largest)
      cycles[i][j]=map<int, unordered_set< vector<int> > >();
      for(int k=(j+1);k<(m.rows()-1);++k){// from the next smallest to the last
        cycles[i][j][k]= unordered_set< vector<int> >();
      }
    }
  }
  int totCycles = 0;
  //cout << m << endl;
  map<int,vector<int> > connected = to_elements(m);
  vector<vector<int> > paths;
  vector<int> conn = connected[input_node-1];
  for(auto it=begin(conn);it!=end(conn);++it){
    vector<int> p;
    p.push_back(input_node-1);
    p.push_back(*it);
    if(*it == (output_node-1)){
      if(save_trajectories)
        trajectories.push_back(p);
    }else{
      paths.push_back(p);
    }
  }
  
  // PARALLEL INFO
  tthread::mutex cycle_mutex;
  tthread::mutex trajectory_mutex;
  tthread::mutex path_mutex;
  NextPathsWorker nextPaths(&cycle_mutex,&trajectory_mutex,&path_mutex,
                            &connected,&paths,&trajectories,&cycles,output_node-1,&aggregated_cycle_participation,
                            &aggregated_trajectory_participation,save_trajectories,&totCycles);
  
  while(paths.size()>0){
    if(use_parallel && paths.size()>limit_for_parallel){
      cout << "WENT PARALLEL : " << paths.size() << endl;
      nextPaths.Restart(&paths);
      parallelFor(0, paths.size(), nextPaths);
      paths = nextPaths.outPaths;
    }else{
      cout << "WENT DIRECT : " << paths.size() << endl;
      paths = next_paths(&connected,&paths,&trajectories,&cycles,output_node-1,&aggregated_cycle_participation,&aggregated_trajectory_participation,save_trajectories,&totCycles);
    }
  }
  vector<vector< vector<int> > > final_out;
  vector< vector<int> > cycles_out;
  for(auto lit=cycles.begin();lit!=cycles.end();++lit){
    for(auto mit=lit->second.begin();mit!=lit->second.end();++mit){
      for(auto sit=mit->second.begin();sit!=mit->second.end();++sit){
        for(auto cycle : sit->second){
          vector<int> c;  
          for(auto n=begin(cycle);n!=end(cycle);n++){
            c.push_back(*n);
          }
          cycles_out.push_back(c);
        }
      }
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

/*Johnson's algorithm is an algorithm for finding all the simple cycles in a directed graph. 
It is based on depth-first search, with modifications in order not to perform checks on the previous 
added cycles each time a new simple cycle is included, and has a worst-case time complexity of O((|V| + |E|)(c + 1)) 
where |V| is the number of vertices in the graph, |E| is the number of edges in the graph, 
and c is the number of cycles in the graph. The main idea is to perform simple cycle searches 
per vertex (in some ordering of the vertices) that always contain that initial vertex so that, 
on the next vertex search, all the edges from and to the previous vertex are removed. 
This reduces any chance of cycle collision on the cumulated list of cycles. 
It also explicits rules to block visited vertices that will not produce any more simple cycles to the 
start vertex during the search, and keep the blockage as long as possible in order to highly reduce 
the amount of ramifications that are to be followed. */

// Function for unblocking vertices
void unblock(int i,vector<bool>& blocked,vector<vector<int> >& B){
  if(blocked[i]){
    blocked[i]=false;
    for(auto it=B[i].begin();it!=B[i].end();++it){
      unblock(*it,blocked,B);
    }
    B[i].clear();
  }
}
// Recursive function to find cycles starting at vertex i 
bool circuit(int starting,int current, int n,vector<int> stack, vector<vector<int> >& A, vector<bool>& blocked, vector<vector<int> >& B,vector<vector<int> >&cycles) 
{
  bool found_cycles = false;
  //Add the current node to the stack and block
  stack.push_back(current);
  blocked[current] = true;
  // Iterate through the ouput edges only from the starting node
  // If the neigbour is blocked skip
  // If it is the same as the first node add simple cycle to output list
  // If it is not blocked continue circuit to that node
  for (int j = starting; j < n; j++){
    if(A[current][j] == 1){
      if(j==starting){ // We've got a simple cycle
        found_cycles = true;
        vector<int> current=stack;
        for(int& x : current) // if you want to add 10 to each element
          x += 1;
        cycles.push_back(current);
        //cycles.push_back(stack);
        if(cycles.size()%10000==0){
          cout << "FAST CYCLES AT : " << cycles.size() << endl;
        }
      }else if(!blocked[j]){
        if(circuit(starting,j,n,stack,A,blocked,B,cycles)){
          found_cycles=true;
        }
      }
    }
  }
  if(found_cycles){
    unblock(current,blocked,B);
  }else{
    for (int j = starting; j < n; j++){
      if(A[current][j] == 1){
        if(find(B[j].begin(),B[j].end(),current)==B[j].end()){
          B[j].push_back(current);
        }
      }
    } 
  }
  return found_cycles;
}
// [[Rcpp::export]]
vector<vector< vector<int> > > fast_cycles(NumericMatrix m){//,bool use_parallel=false){
  NumericMatrix aggregated_cycle_participation( m.rows() );
  vector<vector<int> > cycles; //
  int totCycles = 0;
  vector<vector<int> > A(m.rows(),vector<int>(m.rows()));
  for(int i=0;i<m.rows();++i){
    for(int j=0;j<m.rows();++j){
      A[i][j]=m(i,j);
    }
  }
  // traverse over all vertices 
  for (int i = 0; i < m.rows(); i++){
    // vector to store current cycle 
    vector<int> stack;
    // visited array to keep track of blocked nodes 
    vector<bool>blocked(m.rows(), false);
    // Unblocking array B
    vector<vector<int> > B(m.rows(),vector<int>());
    
    bool dummy = circuit(i,i,m.rows(),stack,A,blocked,B,cycles);
  }
  
  vector<vector< vector<int> > > final_out;
  /*vector<vector<int> > agg_cycles;
  for(int i=0;i<aggregated_cycle_participation.rows();++i){
    vector<int> cycle_row;
    for(int j=0;j<aggregated_cycle_participation.cols();++j){
      cycle_row.push_back(aggregated_cycle_participation(i,j));
    }
    agg_cycles.push_back(cycle_row);
  }*/
  final_out.push_back(cycles); //res[[1]]
  //final_out.push_back(agg_cycles); //res[[2]]
  return final_out;
}

// [[Rcpp::export]]
List fast_cycle_flow(NumericMatrix m,NumericVector x,List cycles){//,bool use_parallel=false){
  NumericMatrix cycle_flow_per_link( m.rows() );
  vector<List> cycle_flows(cycles.size());
  vector<vector<double> > A(m.rows(),vector<double>(m.rows()));
  for(int i=0;i<m.rows();++i){
    for(int j=0;j<m.rows();++j){
      A[i][j]=m(i,j);
      cycle_flow_per_link(i,j)=0.0;
    }
  }
  for(int k=0;k<cycles.size();++k){
    NumericVector cycle = cycles[k];
    int num_nodes = cycle.size();
    double mult = 1.0;
    double stock = 0.0;
    for(int i=0;i<num_nodes-1;++i){
      double link = A[cycle[i]-1][cycle[i+1]-1];
      mult = mult * link;
      stock += x[cycle[i]-1];
    }
    double link = A[cycle[num_nodes-1]-1][cycle[0]-1];
    mult = mult * link;
    stock += x[cycle[num_nodes-1]-1];
    
    double cycle_val = stock*mult;
    double cycle_val_per_link = cycle_val/(double)(cycle.size());
    //cout << "CYCLE : " << k << " NODES : "<< num_nodes << " FLOW : "<< cycle_val << " PER LINK : " << cycle_val_per_link << endl;
    for(int i=0;i<num_nodes-1;++i){
      cycle_flow_per_link(cycle[i]-1,cycle[i+1]-1) = cycle_flow_per_link(cycle[i]-1,cycle[i+1]-1) + cycle_val_per_link;
    }
    cycle_flow_per_link(cycle[num_nodes-1]-1,cycle[0]-1) = cycle_flow_per_link(cycle[num_nodes-1]-1,cycle[0]-1) + cycle_val_per_link;
    cycle_flows[k] = List::create(Named("flow") = cycle_val , _["num_nodes"] = num_nodes);
  }
  List output = List::create(Named("cycle_flow_per_link") = cycle_flow_per_link , _["cycle_flow_per_cycle"] = cycle_flows);
  return output;
}
// FUNCIONES DE CAMBIOS EN CICLOS (recalcular considerando modificaciones en vez de 0)
// ARGUMENTOS DE ENTRADA: 
//   Trayectorias anteriores, Ciclos anteriores, Participacion en Ciclos y Trayectorias anteriores,
//   Links eliminados, Links agregados
// Primer paso es eliminar ciclos y trayectorias que contengan links eliminados. 
//   Recorrer buscando ints seguidos de cada link (entrada y salida del link) en las listas y eliminar.
//   Ojo con ciclos que pueden tener el link al final de la lista (del ultimo elemento al primero)
// Hacer next_paths (nueva funcion de next_paths) del primer link nuevo como paths. En el add_cycles (hacer nueva version) solo agregar ciclos 
// si el ciclo tiene alguno de los links nuevos (sino ya estaba). Checkear si esta el ciclo en la nueva lista igual. Guardar cada link nuevo que se toco durante el recorrido. Guardar las trayectorias desde ese nodo aparte
// La siguiente iteracion se hace solo sobre los links que no se recorrieron durante las busquedas anteriores
// Agregar todos los ciclos nuevos a la lista antigua (si o si tendran links nuevos asi que son nuevos).
// Para las trayectorias buscar para cada link nuevo las trayectorias anteriores que llegaban al nodo de entrada del link
// hacer una lista truncada de estas trayectorias hasta ese nodo. Para cada trayectoria de esa nueva lista, recorrer las nuevas y eliminar las que tocan algun nodo de la antiagua. 
// Las que quedan se concatenan con las nuevas trayectorias y se agregar a las totales. FIN! (ojala jejeje)
/*** R
system.time(res<-paths(matrix(sample(c(0,1),100,replace = TRUE),nrow=10,ncol=10),9,10))
print(length(res[[3]]))
system.time(res<-fast_cycles(matrix(sample(c(0,1),100,replace = TRUE),nrow=10,ncol=10)))
print(length(res[[1]]))
*/
