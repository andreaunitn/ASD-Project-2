#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include "endgame.h"
using namespace std;

//salva in taken l'indice delle pietre inserite nello zaino. ritorna il valore del contenuto dello zaino
int saveKnapSack(int W, vector<int> wt, vector<int> val, int n, vector<bool> &taken);
//salva le pioetre prese con tecnica euristica (rapporto maggiore)
void heuristicKnapSack(int W, vector<pair<double, int> > &rapportoOrdinato, vector<bool> &taken);
//visita il grafo (iviare con true all'inizio)
void tsp(int s, bool partenza, bool mod);
//sceglie se e quali pietre prendere
void takeStones(int s);
//sceglie la prima pietra disponibile
void takeStones2(int s);
//ricostruisce la soluzione
void soluzione();
//stampa la soluzione attuale
void stampa();
//restitutisce se la suluzione attuale è migliore, nel caso la salva
bool bestSol();

struct nodo {
  int pietraPresa = -1; //salva quale pietra è presa in questo nodo, -1 = nessuna presa
  vector<int> adj; //nodi vicini..
  vector<int> dist; //.. e relativa distanza
  vector<int> p; //pietre in questo nodo
  bool visited = false;
};

vector<nodo> grafo;
int N,S;
int M,C;
double R,V_MIN,V_MAX;
vector<int> m; //massa della pietra i
vector<int> e; //energia della pietra i
double E = 0.0, G = 0.0, T = 0.0; //valori da restituire come suluzione
double E1 = -100000000.0, G1 = -100000000.0, T1 = -100000000.0; //soluzione migliore fino ad ora
double peso_tot; //salva il peso totale delle pietre
vector<bool> pietreInZaino; //salva le pietre da inserire nello zaino per indice
vector<int> pietreCitta; //salva dove ogni pietra è stata presa (-1 = non presa)
vector<int> ordineVisite; //salva l'ordine in cui i nodi sono visitati
vector<int> ordineDistanze; //salva le distanze tra nodi nell'ordine in cui sono presi
vector<double> rapporto; //salva il rapporto energia-massa delle pietre
vector<pair<double, int> > rapportoOrdinato; //salva i rapporti ordinati e relativo indice pietra
bool mod = false; //definisce quale takeStones usare -- mod = true -> eseguo takeStones2

ifstream in("input.txt");
ofstream out("output.txt");

////////////////////////////////////////////////////////////////////////////////
//MAIN
int main() {

  in >> N >> S >> M >> C >> R >> V_MIN >> V_MAX;
  grafo.resize(N);
  m.resize(M);
  e.resize(M);
  pietreInZaino.resize(M);
  pietreCitta.resize(M);
  rapporto.resize(M);
  rapportoOrdinato.resize(M);

  if(C >= 100000000 && N > 1950) {
    mod = true;
  }

  for(int i = 0; i < M; i++) {
    int m_i, e_i;
    in >> m_i >> e_i;
    m[i] = m_i;
    e[i] = e_i;
    rapporto[i] = (double)e_i/m_i;
    rapportoOrdinato[i] = make_pair(rapporto[i], i);
    pietreCitta[i] = -1;
    pietreInZaino[i] = mod;
  }

  sort(rapportoOrdinato.begin(), rapportoOrdinato.end());

  for(int i = 0; i < M; i++) {
    int l;
    in >> l;

    for(int j = 0; j < l; j++) {
      int id;
      in >> id;
      grafo[id].p.push_back(i);
    }
  }

  for(int i = 0; i < N; i++) {
    int j = i;
    int t,k=0;
    while(j>0) {
      in >> t;
      grafo[i].dist.push_back(t);
      grafo[i].adj.push_back(k);
      grafo[k].adj.push_back(i);
      grafo[k].dist.push_back(t);
      j--;
      k++;
    }
  }
  //FINE LETTURA INPUT

  if(!mod && C > 50000) {
    heuristicKnapSack(C, rapportoOrdinato, pietreInZaino);
  } else if(!mod){
    saveKnapSack(C, m, e, M, pietreInZaino);
  }

  tsp(S, true, mod);
  soluzione();
  bestSol();
  stampa();

  mod = true;

  heuristicKnapSack(C, rapportoOrdinato, pietreInZaino);

  for(int i = 0; i < M; i++) {
    pietreCitta[i] = -1;
  }

  for(int i = 0; i < N; i++) {
    grafo[i].visited = false;
    grafo[i].pietraPresa = -1;
  }

  peso_tot = 0;
  ordineDistanze.clear();
  ordineVisite.clear();

  E=0.0;
  G=0.0;
  T=0.0;

  tsp(S, true, mod);
  soluzione();
  if(bestSol()) {
    stampa();
  }

  return 0;
}

//FINE MAIN
////////////////////////////////////////////////////////////////////////////////

//funzione accessoria che calcola il tempo impiegato per andare da un nodo
double tempo_impiegato (const int peso, int dist){
  if(C == 0) {return dist/V_MAX;}
  return (double)dist/(V_MAX-((double)peso*(V_MAX-V_MIN)/(double)C));
}

void stampa() {
  out << scientific << setprecision(10) << E1 << " ";
  out << scientific << setprecision(10) << G1 << " ";
  out << scientific << setprecision(10) << T1 << endl;

  for(int i: pietreCitta) {
    out << i << " ";
  }
  out << endl;
  for(int c = ordineVisite.size()-1; c >= 0; c--) {
    out << ordineVisite[c] << " ";
  }
  out << endl;
  out << "***" << endl;
}

void soluzione() {
  double mSoFar = 0.0;
  for(int i = ordineVisite.size()-1; i > 0; i--) {
    if(grafo[ordineVisite[i]].pietraPresa >= 0)
      mSoFar += m[grafo[ordineVisite[i]].pietraPresa];

    T += tempo_impiegato(mSoFar, ordineDistanze[i-1]);
    G += e[grafo[ordineVisite[i]].pietraPresa];
  }
  E = G - R * T;
}

void tsp(int s, bool partenza, bool mod) {
  grafo[s].visited = true;
  if(!partenza) {
    if(mod) {
      takeStones2(s);
    } else {
      takeStones(s);
    }
  }

  ordineVisite.push_back(s);
  unsigned int minSoFar = -1;
  int nodoMin;
  for(int i=0; i<grafo[s].adj.size(); i++) {
    if(!grafo[grafo[s].adj[i]].visited) {
      if(grafo[s].dist[i] < minSoFar) {
        minSoFar = grafo[s].dist[i];
        nodoMin = grafo[s].adj[i];
      }
    }
  }

  if(minSoFar == -1) {
    if(s != S) {
      int res;
      for(int i = 0; i < grafo[S].adj.size(); i++) {
        if(grafo[S].adj[i] == s) {
          res = grafo[S].dist[i];
        }
      }
      ordineDistanze.push_back(res);
      tsp(S, false, mod);
    }
    else {return;}
  } else {
    ordineDistanze.push_back(minSoFar);
    tsp(nodoMin, false, mod);
  }
}

void takeStones(int s) {
  int presa = -1;
  double maxRapp = -1;
  for(int i = 0; i < grafo[s].p.size(); i++) {
    if(pietreInZaino[grafo[s].p[i]]) {
      if(maxRapp < m[grafo[s].p[i]] && m[grafo[s].p[i]] <= C - peso_tot) {
        presa = grafo[s].p[i];
        maxRapp = rapporto[grafo[s].p[i]];
      }
    }
  }
  if(presa == -1) return;
  peso_tot += m[presa];
  pietreInZaino[presa] = false;
  pietreCitta[presa] = s;
  grafo[s].pietraPresa = presa;
}

void takeStones2(int s) {
  int presa = -1;
  for(int i = 0; i < grafo[s].p.size(); i++) {
    if(pietreInZaino[grafo[s].p[i]]) {
      if(m[grafo[s].p[i]] <= C - peso_tot) {
        presa = grafo[s].p[i];
      }
    }
  }

  if(presa == -1) return;
  peso_tot += m[presa];
  pietreInZaino[presa] = false;
  pietreCitta[presa] = s;
  grafo[s].pietraPresa = presa;
}

int saveKnapSack(int W, vector<int> wt, vector<int> val, int n, vector<bool> &taken) {

  int i, w;
  int maxPr;
  int **K = new int*[n+1];
  for(int i = 0; i < n+1; i++) {
    K[i] = new int[W+1];
  }
  for (i = 0; i <= n; i++) {
    for (w = 0; w <= W; w++) {
      if (i == 0 || w == 0)
        K[i][w] = 0;
      else if (wt[i - 1] <= w)
        K[i][w] = max(val[i - 1] + K[i - 1][w - wt[i - 1]], K[i - 1][w]);
      else
        K[i][w] = K[i - 1][w];
      }
  }

  int res = K[n][W];
  maxPr = res;

  w = W;
  for (i = n; i > 0 && res > 0; i--) {
    if (res == K[i - 1][w])
      continue;
    else {
      taken[i - 1] = true;
      res = res - val[i - 1];
      w = w - wt[i - 1];
    }
  }

  for(int i = 0; i < n+1; i++) {
    delete []K[i];
  }

  delete []K;

  return maxPr;
}

void heuristicKnapSack(int W, vector<pair<double, int> > &rapportoOrdinato, vector<bool> &taken) {
  int w = 0;
  for(int i=rapportoOrdinato.size()-1; i>0; i--) {
    if(m[rapportoOrdinato[i].second] <= W - w) {
      taken[rapportoOrdinato[i].second] = true;
      w += m[rapportoOrdinato[i].second];
    }
  }
}

bool bestSol() {
  bool res = false;
  if(E > E1) {
    E1 = E;
    G1 = G;
    T1 = T;
    res = true;
  }
  return res;
}
