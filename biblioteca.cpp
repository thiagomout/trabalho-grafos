// g++ -O2 -std=c++17 biblioteca.cpp -o estudos
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <tuple>
#include <algorithm>
#include <queue>
#include <stack>
#include <map>
#include <iomanip>
#include <climits>
#include <memory>
#include <chrono>
#include <numeric> 

using namespace std;

// -------------------- Interface Abstrata do Grafo --------------------
class Grafo {
public:
    virtual ~Grafo() = default;
    virtual int get_num_vertices() const = 0;
    virtual long long get_num_arestas() const = 0;
    virtual vector<int> get_vizinhos(int u) const = 0;
};

// -------------------- Implementação com Lista de Adjacência --------------------
class GrafoListaAdj : public Grafo {
private:
    int n;
    long long m;
    vector<vector<int>> adj;

public:
    GrafoListaAdj(int num_vertices, long long num_arestas, vector<vector<int>>&& lista_adj)
        : n(num_vertices), m(num_arestas), adj(move(lista_adj)) {}

    int get_num_vertices() const override { return n; }
    long long get_num_arestas() const override { return m; }
    vector<int> get_vizinhos(int u) const override {
        if (u > 0 && (size_t)u < adj.size()) return adj[u];
        return {};
    }
};

// -------------------- Implementação com Matriz de Adjacência --------------------
class GrafoMatrizAdj : public Grafo {
private:
    int n;
    long long m;
    vector<vector<int>> M;

public:
    GrafoMatrizAdj(int num_vertices, long long num_arestas, const vector<vector<int>>& adj_list)
        : n(num_vertices), m(num_arestas) {
        M.assign(n, vector<int>(n, 0));
        for (int u = 1; u <= n; ++u) {
            if((size_t)u < adj_list.size()) {
                for (int v : adj_list[u]) {
                    M[u - 1][v - 1] = 1;
                }
            }
        }
    }

    int get_num_vertices() const override { return n; }
    long long get_num_arestas() const override { return m; }
    vector<int> get_vizinhos(int u) const override {
        vector<int> vizinhos;
        if (u > 0 && u <= n) {
            for (int j = 0; j < n; ++j) {
                if (M[u - 1][j] == 1) {
                    vizinhos.push_back(j + 1);
                }
            }
        }
        return vizinhos;
    }
};

// -------------------- Leitura do grafo --------------------
tuple<int, long long, vector<vector<int>>> ler_dados_grafo(const string& caminho) {
    ifstream f(caminho);
    if (!f) {
        cerr << "Erro: Nao foi possivel abrir o arquivo: " << caminho << "\n";
        exit(1);
    }
    int n;
    f >> n;
    string rest;
    getline(f, rest);

    vector<vector<int>> adj(n + 1);
    long long m = 0;
    int u, v;
    while (f >> u >> v) {
        if (u > n || v > n || u < 1 || v < 1) continue;
        adj[u].push_back(v);
        adj[v].push_back(u);
        m++;
    }
    return {n, m, adj};
}

// -------------------- Algoritmos Genéricos --------------------

pair<vector<int>, vector<int>> bfs(const Grafo& grafo, int s) {
    int n = grafo.get_num_vertices();
    vector<int> pai(n + 1, -1), nivel(n + 1, -1);
    if (s < 1 || s > n) return {pai, nivel};
    
    queue<int> fila;
    nivel[s] = 0;
    fila.push(s);

    while (!fila.empty()) {
        int u = fila.front();
        fila.pop();
        for (int v : grafo.get_vizinhos(u)) {
            if (nivel[v] == -1) {
                nivel[v] = nivel[u] + 1;
                pai[v] = u;
                fila.push(v);
            }
        }
    }
    return {pai, nivel};
}

pair<vector<int>, vector<int>> dfs(const Grafo& grafo, int s) {
    int n = grafo.get_num_vertices();
    vector<int> pai(n + 1, -1), nivel(n + 1, -1);
    if (s < 1 || s > n) return {pai, nivel};

    stack<pair<int, int>> pilha; // {vértice, nível}
    vector<bool> visitado(n + 1, false);

    pilha.push({s, 0});
    visitado[s] = true;
    nivel[s] = 0;

    while (!pilha.empty()) {
        auto [u, niv_u] = pilha.top();
        pilha.pop();
        
        auto vizinhos = grafo.get_vizinhos(u);
        reverse(vizinhos.begin(), vizinhos.end());

        for (int v : vizinhos) {
            if (!visitado[v]) {
                visitado[v] = true;
                pai[v] = u;
                nivel[v] = niv_u + 1;
                pilha.push({v, niv_u + 1});
            }
        }
    }
    return {pai, nivel};
}

int distancia(const Grafo& grafo, int u, int v) {
    auto [_, nivel] = bfs(grafo, u);
    return nivel[v];
}

vector<vector<int>> componentes_conexas(const Grafo& grafo) {
    int n = grafo.get_num_vertices();
    if (n == 0) return {};
    vector<bool> visitado(n + 1, false);
    vector<vector<int>> comps;

    for (int i = 1; i <= n; ++i) {
        if (!visitado[i]) {
            vector<int> comp_atual;
            queue<int> q;
            q.push(i);
            visitado[i] = true;

            while(!q.empty()){
                int u = q.front();
                q.pop();
                comp_atual.push_back(u);
                for(int v : grafo.get_vizinhos(u)) {
                    if(!visitado[v]){
                        visitado[v] = true;
                        q.push(v);
                    }
                }
            }
            comps.push_back(move(comp_atual));
        }
    }

    sort(comps.begin(), comps.end(), [](const auto& a, const auto& b) {
        return a.size() > b.size();
    });
    return comps;
}

int diametro_exato(const Grafo& grafo) {
    int n = grafo.get_num_vertices();
    int diam = 0;
    for (int i = 1; i <= n; ++i) {
        auto [_, nivel] = bfs(grafo, i);
        for (int j = 1; j <= n; ++j) {
            if (nivel[j] == -1) return INT_MAX; 
            if (nivel[j] > diam) diam = nivel[j];
        }
    }
    return diam;
}

int diametro_aproximado(const Grafo& grafo) {
    int n = grafo.get_num_vertices();
    if (n <= 1) return 0;
    
    int a = 1;
    auto [p1, n1] = bfs(grafo, a);
    int b = a;
    for(int i = 1; i <= n; ++i) {
        if (n1[i] > n1[b]) b = i;
    }

    auto [p2, n2] = bfs(grafo, b);
    int diam_aprox = 0;
    for(int i = 1; i <= n; ++i) {
        if (n2[i] > diam_aprox) diam_aprox = n2[i];
    }
    return diam_aprox;
}

tuple<int, int, double, double> estatisticas_grau(const Grafo& grafo) {
    int n = grafo.get_num_vertices();
    if (n == 0) return {0, 0, 0.0, 0.0};
    
    vector<int> graus(n);
    long long soma_graus = 0;
    for (int i = 0; i < n; ++i) {
        int grau = grafo.get_vizinhos(i + 1).size();
        graus[i] = grau;
        soma_graus += grau;
    }
    
    sort(graus.begin(), graus.end());
    
    int gmin = graus.front();
    int gmax = graus.back();
    double media = (double)soma_graus / n;
    double mediana = (n % 2 == 0) ? (graus[n/2 - 1] + graus[n/2]) / 2.0 : graus[n/2];
    
    return {gmin, gmax, media, mediana};
}

void pausa(const string& msg) {
    cout << endl << msg << endl;
    cout << "Pressione ENTER para continuar...";
    string dummy;
    getline(cin, dummy);
    cout << endl;
}

// -------------------- Driver --------------------
int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc != 2) {
        cerr << "Erro: Uso incorreto.\nExecute como: ./estudos <lista|matriz>\n";
        return 1;
    }
    string tipo_representacao = argv[1];
    if (tipo_representacao != "lista" && tipo_representacao != "matriz") {
        cerr << "Erro: Representacao '" << tipo_representacao << "' invalida. Use 'lista' ou 'matriz'.\n";
        return 1;
    }
    
    const string caminho_arquivo = "grafo_2.txt";
    const string arquivo_saida = "RelatorioCompleto.txt";

    cout << "Lendo dados do grafo de '" << caminho_arquivo << "'...\n";
    auto [n, m, adj_data] = ler_dados_grafo(caminho_arquivo);
    
    unique_ptr<Grafo> grafo;
    if (tipo_representacao == "lista") {
        grafo = make_unique<GrafoListaAdj>(n, m, move(adj_data));
        cout << "Grafo carregado com representacao em LISTA de adjacencia.\n";
    } else {
        grafo = make_unique<GrafoMatrizAdj>(n, m, adj_data);
        cout << "Grafo carregado com representacao em MATRIZ de adjacencia.\n";
    }

    pausa(">>> O grafo foi carregado na memoria. Verifique o consumo do processo agora. <<<");
    
    cout << "\nIniciando benchmark de performance das buscas...\n";
    const int N_RUNS = 100;
    vector<int> seeds;
    if (n > 0) {
        for(int i = 0; i < N_RUNS; ++i) {
            seeds.push_back((i * (n / N_RUNS)) % n + 1);
        }
    }
    
    auto bfs_start = chrono::high_resolution_clock::now();
    for (int s : seeds) {
        volatile auto r = bfs(*grafo, s);
    }
    auto bfs_end = chrono::high_resolution_clock::now();
    auto bfs_duration = chrono::duration_cast<chrono::milliseconds>(bfs_end - bfs_start);
    double bfs_avg = (double)bfs_duration.count() / N_RUNS;
    
    auto dfs_start = chrono::high_resolution_clock::now();
    for (int s : seeds) {
        volatile auto r = dfs(*grafo, s);
    }
    auto dfs_end = chrono::high_resolution_clock::now();
    auto dfs_duration = chrono::duration_cast<chrono::milliseconds>(dfs_end - dfs_start);
    double dfs_avg = (double)dfs_duration.count() / N_RUNS;
    
    cout << "Benchmark concluido.\n";

    cout << "\nIniciando analise completa do grafo...\n";
    
    auto [gmin, gmax, gmed, gmediana] = estatisticas_grau(*grafo);
    auto comps = componentes_conexas(*grafo);
    int diam_aprox = diametro_aproximado(*grafo);

    int diam_exato = -1;
    const int LIMITE_N_DIAMETRO = 2000;
    if (n > 0 && n <= LIMITE_N_DIAMETRO) {
        cout << "Calculando diametro exato (n=" << n << ", pode demorar um pouco)...\n";
        diam_exato = diametro_exato(*grafo);
    } else {
        cout << "AVISO: Grafo muito grande (n=" << n << "). Pulando calculo do diametro exato.\n";
    }
    
    vector<int> origens_pais = {1, 2, 3};
    vector<int> alvos_pais = {10, 20, 30};
    map<int, vector<int>> bfs_pais_resultados, dfs_pais_resultados;
    for (int s : origens_pais) {
        auto [pais_b, _b] = bfs(*grafo, s);
        auto [pais_d, _d] = dfs(*grafo, s);
        for (int alvo : alvos_pais) {
            bfs_pais_resultados[s].push_back(pais_b[alvo]);
            dfs_pais_resultados[s].push_back(pais_d[alvo]);
        }
    }
    
    vector<pair<int, int>> pares_dist = {{10, 20}, {10, 30}, {20, 30}};
    vector<int> resultados_dist;
    for (const auto& par : pares_dist) {
        resultados_dist.push_back(distancia(*grafo, par.first, par.second));
    }

    cout << "Analise concluida. Gerando relatorio...\n";

    ofstream out(arquivo_saida);
    out << "=== RELATORIO DE ANALISE DE GRAFO ===\n\n";
    out << "Arquivo de entrada: " << caminho_arquivo << "\n";
    out << "Representacao utilizada: " << tipo_representacao << "\n\n";
    
    out << "--- Informacoes Basicas ---\n";
    out << "Vertices: " << n << "\n";
    out << "Arestas: " << m << "\n\n";
    
    out << "--- Benchmark de Performance (" << N_RUNS << " execucoes) ---\n";
    out << "Tempo total BFS: " << bfs_duration.count() << " ms\n";
    out << "Tempo medio por BFS: " << fixed << setprecision(4) << bfs_avg << " ms\n";
    out << "Tempo total DFS: " << dfs_duration.count() << " ms\n";
    out << "Tempo medio por DFS: " << fixed << setprecision(4) << dfs_avg << " ms\n\n";

    out << "--- Estatisticas de Grau ---\n";
    out << "Grau Minimo: " << gmin << "\n";
    out << "Grau Maximo: " << gmax << "\n";
    out << "Grau Medio: " << gmed << "\n";
    out << "Mediana do Grau: " << gmediana << "\n\n";
    
    out << "--- Componentes Conexas ---\n";
    out << "Quantidade: " << comps.size() << "\n";
    out << "Tamanho da Maior Componente: " << (comps.empty() ? 0 : comps.front().size()) << "\n";
    out << "Tamanho da Menor Componente: " << (comps.empty() ? 0 : comps.back().size()) << "\n\n";
    
    out << "--- Diametro ---\n";
    if (diam_exato != -1) {
        out << "Exato: " << (diam_exato == INT_MAX ? "Infinito (grafo desconexo)" : to_string(diam_exato)) << "\n";
    } else {
        out << "Exato: Nao calculado (grafo muito grande)\n";
    }
    out << "Aproximado (2-BFS): " << diam_aprox << "\n\n";
    
    out << "--- Item 4.4: Pais dos Vertices 10, 20, 30 ---\n";
    for (int s : origens_pais) {
        out << "Buscas a partir do vertice " << s << ":\n";
        out << "  [BFS] Pais de (10, 20, 30): (" << bfs_pais_resultados[s][0] << ", " << bfs_pais_resultados[s][1] << ", " << bfs_pais_resultados[s][2] << ")\n";
        out << "  [DFS] Pais de (10, 20, 30): (" << dfs_pais_resultados[s][0] << ", " << dfs_pais_resultados[s][1] << ", " << dfs_pais_resultados[s][2] << ")\n";
    }
    out << "\n";

    out << "--- Item 4.5: Distancia entre Pares ---\n";
    out << "  Distancia (10, 20): " << resultados_dist[0] << "\n";
    out << "  Distancia (10, 30): " << resultados_dist[1] << "\n";
    out << "  Distancia (20, 30): " << resultados_dist[2] << "\n";

    out.close();
    cout << "\nRelatorio salvo em '" << arquivo_saida << "'.\n";

    return 0;
}