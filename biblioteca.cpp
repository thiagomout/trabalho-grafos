// g++ -O2 -std=c++17 biblioteca.cpp -o biblioteca
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
#include <limits>
#include <queue>
#include <random>
#include <filesystem> 

using namespace std;

// -------------------- Interface Abstrata do Grafo --------------------
class Grafo {
public:
    virtual ~Grafo() = default;
    virtual int get_num_vertices() const = 0;
    virtual long long get_num_arestas() const = 0;
    virtual vector<int> get_vizinhos(int u) const = 0;
    virtual vector<pair<int, double>> get_vizinhos_com_peso(int u) const = 0;
};

// -------------------- Implementação com Lista de Adjacência --------------------
class GrafoListaAdj : public Grafo {
private:
    int n;
    long long m;
    vector<vector<pair<int, double>>> adj;

public:
    GrafoListaAdj(int num_vertices, long long num_arestas, vector<vector<pair<int, double>>>&& lista_adj)
        : n(num_vertices), m(num_arestas), adj(move(lista_adj)) {}

    int get_num_vertices() const override { return n; }
    long long get_num_arestas() const override { return m; }
    vector<int> get_vizinhos(int u) const override {
        vector<int> vizinhos_ids;
        if (u > 0 && (size_t)u < adj.size()) {
            vizinhos_ids.reserve(adj[u].size());
            for (const auto& par : adj[u]) {
                vizinhos_ids.push_back(par.first);
            }
        }
        return vizinhos_ids;
    }

    vector<pair<int, double>> get_vizinhos_com_peso(int u) const override {
        if (u > 0 && (size_t)u < adj.size()) {
            return adj[u];
        }
        return {};
    }
};

// -------------------- Implementação com Matriz de Adjacência --------------------
class GrafoMatrizAdj : public Grafo {
private:
    int n;
    long long m;
    vector<vector<double>> M;

public:
    GrafoMatrizAdj(int num_vertices, long long num_arestas, const vector<vector<pair<int, double>>>& adj_list)
        : n(num_vertices), m(num_arestas) {
        const double INF = std::numeric_limits<double>::infinity();
        M.assign(n, vector<double>(n, INF));
        for (int u = 1; u <= n; ++u) {
            if((size_t)u < adj_list.size()) {
                for (const auto& par : adj_list[u]) {
                    int v = par.first;
                    double w = par.second;
                    M[u - 1][v - 1] = w;
                }
            }
        }
    }

    int get_num_vertices() const override { return n; }
    long long get_num_arestas() const override { return m; }
    vector<int> get_vizinhos(int u) const override {
        vector<int> vizinhos_ids;
        if (u > 0 && u <= n) {
            const double INF = std::numeric_limits<double>::infinity();
            for (int j = 0; j < n; ++j) {
                if (M[u - 1][j] != INF) {
                    vizinhos_ids.push_back(j + 1);
                }
            }
        }
        return vizinhos_ids;
    }
    vector<pair<int, double>> get_vizinhos_com_peso(int u) const override {
        vector<pair<int, double>> vizinhos_com_peso;
        if (u > 0 && u <= n) {
            const double INF = std::numeric_limits<double>::infinity();
            for (int j = 0; j < n; ++j) {
                if (M[u - 1][j] != INF) {
                    vizinhos_com_peso.push_back({j + 1, M[u - 1][j]});
                }
            }
        }
        return vizinhos_com_peso;
    }
};


// -------------------- Leitura de dados --------------------

string trim(const string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (string::npos == first) {
        return str;
    }
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, (last - first + 1));
}

// [MODIFICADO] Adicionado parâmetro 'direcionado'
tuple<int, long long, vector<vector<pair<int, double>>>, bool> ler_dados_grafo(const string& caminho, bool direcionado) {
    ifstream f(caminho);
    if (!f) {
        cerr << "Erro: Nao foi possivel abrir o arquivo: " << caminho << "\n";
        exit(1);
    }
    int n;
    f >> n;
    string rest;
    getline(f, rest);

    vector<vector<pair<int, double>>> adj(n + 1);
    long long m = 0;
    int u, v;
    double w;

    bool tem_peso_negativo = false;

    while (f >> u >> v >> w) {
        if (u > n || v > n || u < 1 || v < 1) continue;

        if (w < 0) tem_peso_negativo = true;
        
        // Aresta de ida (u -> v) sempre existe
        adj[u].push_back({v, w});
        
        // [MODIFICADO] Só adiciona a volta se NÃO for direcionado
        if (!direcionado) {
            adj[v].push_back({u, w});
        }
        
        m++; // Conta a aresta/arco
    }
    return {n, m, adj, tem_peso_negativo};
}

// [RE-MODIFICADO] Versão mais robusta para ler 'ID,NOME'
pair<map<string, int>, vector<string>> ler_mapeamento_vertices(const string& caminho_mapa, int n) {
    ifstream f(caminho_mapa);
    if (!f) {
        cerr << "Erro: Nao foi possivel abrir o arquivo de mapeamento: " << caminho_mapa << "\n";
        return {{}, vector<string>(n + 1, "N/A")};
    }

    map<string, int> nome_para_id;
    vector<string> id_para_nome(n + 1, "NOME_DESCONHECIDO"); 
    string linha;

    // cout << "\n--- DEBUG: INICIANDO LEITURA (Modo Robusto) ---" << endl;

    char bom[3];
    f.read(bom, 3);
    if (f.gcount() == 3 && (unsigned char)bom[0] == 0xEF && (unsigned char)bom[1] == 0xBB && (unsigned char)bom[2] == 0xBF) {
        // BOM encontrado e ignorado
    } else {
        f.seekg(0);
    }

    int linha_count = 0;
    while (getline(f, linha)) {
        linha_count++;
        if (linha.empty() || linha == "\r") continue; 

        size_t pos_virgula = linha.find(',');
        if (pos_virgula == string::npos) {
            continue;
        }

        string id_str = linha.substr(0, pos_virgula);
        string nome_str = linha.substr(pos_virgula + 1);

        int id;
        try {
            id = stoi(trim(id_str));
        } catch (...) {
            continue;
        }

        string nome_final = trim(nome_str);
        
        if (id > 0 && id <= n) {
            nome_para_id[nome_final] = id;
            id_para_nome[id] = nome_final;
        }
    }
    
    cout << "Mapeamento lido. " << nome_para_id.size() << " nomes carregados de " << linha_count << " linhas.\n";
    return {nome_para_id, id_para_nome};
}

// -------------------- Algoritmos Genéricos --------------------

class FilaDePrioridade {
public:
    virtual ~FilaDePrioridade() = default;
    virtual void inicializar(int num_vertices, int s) = 0;
    virtual bool esta_vazia() const = 0;
    virtual int extrair_min() = 0;
    virtual void inserir_ou_atualizar(int vertice, double nova_distancia) = 0;
    virtual double get_distancia(int v) const = 0;
};

class FilaComVetor : public FilaDePrioridade {
private:
    vector<double> dist;
    vector<bool> visitado;
    int n;
    int restantes; 

public:
    void inicializar(int num_vertices, int s) override {
        n = num_vertices;
        restantes = n;
        const double INFINITO = std::numeric_limits<double>::infinity();
        dist.assign(n + 1, INFINITO);
        visitado.assign(n + 1, false);
        dist[s] = 0.0;
    }

    bool esta_vazia() const override {
        return restantes == 0;
    }

    int extrair_min() override {
        double min_dist = std::numeric_limits<double>::infinity();
        int min_v = -1;
        for (int v = 1; v <= n; ++v) {
            if (!visitado[v] && dist[v] < min_dist) {
                min_dist = dist[v];
                min_v = v;
            }
        }
        if (min_v != -1) {
            visitado[min_v] = true; 
            restantes--;
        }
        return min_v; 
    }

    void inserir_ou_atualizar(int vertice, double nova_distancia) override {
        dist[vertice] = nova_distancia;
    }
    double get_distancia(int v) const override {
        return dist[v];
    }
};

class FilaComHeap : public FilaDePrioridade {
private:
    using Par = pair<double, int>; 
    std::priority_queue<Par, vector<Par>, greater<Par>> heap;
    vector<double> dist;
    int n;
public:
    void inicializar(int num_vertices, int s) override {
        n = num_vertices;
        const double INFINITO = std::numeric_limits<double>::infinity();
        dist.assign(n + 1, INFINITO);
        dist[s] = 0.0;
        heap.push({0.0, s});
    }

    bool esta_vazia() const override {
        return heap.empty();
    }

    int extrair_min() override {
        while (!heap.empty()) {
            auto [d, v] = heap.top();
            heap.pop();
            if (d == dist[v]) {
                return v;
            }
        }
        return -1; 
    }
    void inserir_ou_atualizar(int vertice, double nova_distancia) override {
        if (nova_distancia < dist[vertice]) {
            dist[vertice] = nova_distancia;
            heap.push({nova_distancia, vertice});
        }
    }

    double get_distancia(int v) const override {
        return dist[v];
    }
};

pair<vector<double>, vector<int>> dijkstra(const Grafo& grafo, int s, FilaDePrioridade& fila) {
    int n = grafo.get_num_vertices();
    const double INFINITO = std::numeric_limits<double>::infinity();
    vector<int> pai(n + 1, -1);

    fila.inicializar(n, s);

    while (!fila.esta_vazia()) {
        int u = fila.extrair_min();
        if (u == -1 || fila.get_distancia(u) == INFINITO) break;

        for (const auto& par : grafo.get_vizinhos_com_peso(u)) {
            int v = par.first;
            double peso = par.second;
            double nova_dist = fila.get_distancia(u) + peso;
            if (nova_dist < fila.get_distancia(v)) {
                fila.inserir_ou_atualizar(v, nova_dist);
                pai[v] = u;
            }
        }
    }

    vector<double> dist_final(n + 1, INFINITO);
    for (int v = 1; v <= n; ++v) {
        dist_final[v] = fila.get_distancia(v);
    }
    return {dist_final, pai};
}

/**
 * @brief Bellman-Ford (SPFA - Shortest Path Faster Algorithm).
 * Calcula distâncias de TODOS os vértices para um DESTINO 't'.
 * Implementa as Otimizações 1 e 2 solicitadas.
 */
pair<bool, pair<vector<double>, vector<int>>> bellman_ford(const Grafo& grafo, int t) {
    int n = grafo.get_num_vertices();
    const double INFINITO = std::numeric_limits<double>::infinity();
    
    // Estruturas de retorno
    vector<double> dist(n + 1, INFINITO);
    vector<int> pai(n + 1, -1); 
    
    // 1. Construção do Grafo Invertido
    // Necessário porque queremos ir de TODOS para UM. O algoritmo roda "de trás para frente".
    vector<vector<pair<int, double>>> adj_inv(n + 1);
    for (int u = 1; u <= n; ++u) {
        for (auto& par : grafo.get_vizinhos_com_peso(u)) {
            int v = par.first;
            double w = par.second;
            // Se existe aresta u->v original, guardamos v->u no inverso
            adj_inv[v].push_back({u, w});
        }
    }

    // 2. Inicialização das Estruturas do SPFA
    queue<int> fila;
    vector<bool> na_fila(n + 1, false); 
    
    // Vetor para controlar a profundidade do caminho (para detectar ciclo negativo)
    // Se um caminho tiver mais que N arestas, existe um ciclo.
    vector<int> arestas_no_caminho(n + 1, 0); 

    // Inicializa o destino (que age como fonte na busca reversa)
    dist[t] = 0.0;
    fila.push(t);
    na_fila[t] = true;

    // 3. Loop Principal (SPFA)
    while (!fila.empty()) {
        int v = fila.front(); 
        fila.pop();
        na_fila[v] = false;

        // Relaxamento Reverso: Olhamos quem chega em 'v' (u -> v)
        for (auto& par : adj_inv[v]) {
            int u = par.first;     // Vizinho que aponta para v
            double peso = par.second; // Peso da aresta u->v

            // Se achamos um caminho melhor para 't' saindo de 'u' passando por 'v'
            if (dist[u] > peso + dist[v]) {
                dist[u] = peso + dist[v];
                pai[u] = v; // O "pai" de u é v (pois estamos indo em direção a t)
                
                // Detecção de Ciclo Negativo:
                // O caminho até u agora tem 1 aresta a mais que o caminho até v.
                arestas_no_caminho[u] = arestas_no_caminho[v] + 1;
                
                if (arestas_no_caminho[u] > n) {
                     return {true, {dist, pai}}; // Ciclo negativo detectado
                }

                // Se u não está na fila para ser processado, adiciona
                if (!na_fila[u]) {
                    fila.push(u);
                    na_fila[u] = true;
                }
            }
        }
    }

    return {false, {dist, pai}};
}

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

    stack<pair<int, int>> pilha; 
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

namespace RelatorioUtils {
    string get_caminho_str(int origem, int destino, const vector<int>& pai) {
        if (pai[destino] == -1 && destino != origem) {
            return (destino == origem) ? to_string(origem) : "Nenhum";
        }
        stringstream ss;
        stack<int> caminho;
        int atual = destino;
        while (atual != -1) { 
            caminho.push(atual);
            atual = pai[atual];
        }
        while (!caminho.empty()) {
            ss << caminho.top();
            caminho.pop();
            if (!caminho.empty()) {
                ss << " -> ";
            }
        }
        return ss.str();
    }

    string get_caminho_str_com_nomes(int origem_id, int destino_id, const vector<int>& pai, const vector<string>& id_para_nome) {
        if (pai[destino_id] == -1 && destino_id != origem_id) {
            return "Nenhum";
        }
        stringstream ss;
        stack<int> caminho_ids;
        int atual = destino_id;
        while (atual != -1) {
            caminho_ids.push(atual);
            atual = pai[atual];
        }
        while (!caminho_ids.empty()) {
            int id_atual = caminho_ids.top();
            caminho_ids.pop();
            string nome = (id_atual > 0 && (size_t)id_atual < id_para_nome.size()) 
                            ? id_para_nome[id_atual] 
                            : ("ID(" + to_string(id_atual) + ")");
            ss << nome;
            if (!caminho_ids.empty()) {
                ss << " -> ";
            }
        }
        return ss.str();
    }

    void imprimir_linha_tabela(ostream& out, const string& col1, const string& col2, const string& col3) {
        out << "   | " << left << setw(10) << col1
            << " | " << left << setw(18) << col2
            << " | " << left << setw(35) << col3 << " |" << endl;
    }

    void imprimir_borda_tabela_dijkstra(ostream& out) {
        out << "   +------------+--------------------+-------------------------------------+" << endl;
    }

    void imprimir_linha_tabela_nomes(ostream& out, const string& col1, const string& col2, const string& col3) {
        out << " | " << left << setw(25) << col1 
            << " | " << left << setw(18) << col2 
            << " | " << left << setw(50) << col3 << " |" << endl; 
    }

    void imprimir_borda_tabela_nomes(ostream& out) {
        out << " +---------------------------+--------------------+----------------------------------------------------+" << endl;
    }

    void imprimir_tabela_benchmark(ostream& out, 
                                     long long t_total_vetor_ms, double t_medio_vetor, 
                                     long long t_total_heap_ms, double t_medio_heap) 
    {
        out << "   +------------------+-------------------+------------------+" << endl;
        out << "   | Implementacao    | Tempo Total (ms)  | Tempo Medio (ms) |" << endl;
        out << "   +------------------+-------------------+------------------+" << endl;
        out << fixed << setprecision(4); 
        
        if (t_total_vetor_ms > 0 || t_medio_vetor > 0.0) {
            out << "   | Vetor (O(V^2))   | " << left << setw(17) << t_total_vetor_ms 
                << " | " << left << setw(16) << t_medio_vetor << " |" << endl;
        }
        if (t_total_heap_ms > 0 || t_medio_heap > 0.0) {
            out << "   | Heap (O(E log V))| " << left << setw(17) << t_total_heap_ms 
                << " | " << left << setw(16) << t_medio_heap << " |" << endl;
        }
        out << "   +------------------+-------------------+------------------+" << endl;
    }
}

// -------------------- Driver --------------------
int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc < 4) {
        cerr << "Erro: Uso incorreto.\n";
        cerr << "Modo Relatorio: ./estudos <lista|matriz> <arquivo.txt> relatorio [direcionado]\n";
        cerr << "Modo Dijkstra:  ./estudos <lista|matriz> <arquivo.txt> dijkstra [k_benchmark] [vetor|heap|ambos] [direcionado]\n";
        return 1;
    }
    string tipo_representacao = argv[1];
    if (tipo_representacao != "lista" && tipo_representacao != "matriz") {
        cerr << "Erro: Representacao '" << tipo_representacao << "' invalida. Use 'lista' ou 'matriz'.\n";
        return 1;
    }
    string caminho_arquivo = argv[2];
    string modo = argv[3];

    // --- [NOVO] Verifica se é direcionado (Tarefa 1) ---
    bool eh_direcionado = false;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "direcionado") {
            eh_direcionado = true;
            break;
        }
    }

    string nome_base_arquivo;
    size_t pos_barra = caminho_arquivo.find_last_of("/\\");
    if (pos_barra != string::npos) {
        nome_base_arquivo = caminho_arquivo.substr(pos_barra + 1);
    } else {
        nome_base_arquivo = caminho_arquivo;
    }
    size_t pos_ponto = nome_base_arquivo.find_last_of(".");
    if (pos_ponto != string::npos) {
        nome_base_arquivo = nome_base_arquivo.substr(0, pos_ponto);
    }


    // --- 2. Leitura do Grafo ---
    cout << "Lendo dados do grafo de '" << caminho_arquivo << "'...\n";
    if (eh_direcionado) {
        cout << "-> MODO DIRECIONADO ATIVADO.\n";
    } else {
        cout << "-> MODO NAO-DIRECIONADO (Padrao).\n";
    }

    // [MODIFICADO] Passa o flag eh_direcionado
    auto [n, m, adj_data, tem_peso_negativo] = ler_dados_grafo(caminho_arquivo, eh_direcionado);
    
    unique_ptr<Grafo> grafo;
    if (tipo_representacao == "lista") {
        grafo = make_unique<GrafoListaAdj>(n, m, move(adj_data));
        cout << "Grafo carregado com representacao em LISTA de adjacencia.\n";
    } else {
        grafo = make_unique<GrafoMatrizAdj>(n, m, adj_data);
        cout << "Grafo carregado com representacao em MATRIZ de adjacencia.\n";
    }
    cout << "Vertices: " << n << ", Arestas/Arcos: " << m << endl;


    if (modo == "relatorio") {
        const string arquivo_saida = "Relatorio_" + nome_base_arquivo + ".txt";
        
        pausa(">>> O grafo foi carregado na memoria. Verifique o consumo do processo agora. <<<");
    
        cout << "\nIniciando benchmark de performance das buscas (BFS/DFS)...\n";
        const int N_RUNS = 100;
        vector<int> seeds;
        if (n > 0) {
            for(int i = 0; i < N_RUNS; ++i) {
                seeds.push_back((i * (n / max(1, N_RUNS))) % n + 1);
            }
        }
        
        auto bfs_start = chrono::high_resolution_clock::now();
        for (int s : seeds) { volatile auto r = bfs(*grafo, s); }
        auto bfs_end = chrono::high_resolution_clock::now();
        auto bfs_duration = chrono::duration_cast<chrono::milliseconds>(bfs_end - bfs_start);
        double bfs_avg = (N_RUNS > 0) ? (double)bfs_duration.count() / N_RUNS : 0.0;
        
        auto dfs_start = chrono::high_resolution_clock::now();
        for (int s : seeds) { volatile auto r = dfs(*grafo, s); }
        auto dfs_end = chrono::high_resolution_clock::now();
        auto dfs_duration = chrono::duration_cast<chrono::milliseconds>(dfs_end - dfs_start);
        double dfs_avg = (N_RUNS > 0) ? (double)dfs_duration.count() / N_RUNS : 0.0;
        
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
            if (s > n) continue;
            auto [pais_b, _b] = bfs(*grafo, s);
            auto [pais_d, _d] = dfs(*grafo, s);
            for (int alvo : alvos_pais) {
                if(alvo > n) continue;
                bfs_pais_resultados[s].push_back(pais_b[alvo]);
                dfs_pais_resultados[s].push_back(pais_d[alvo]);
            }
        }
        
        vector<pair<int, int>> pares_dist = {{10, 20}, {10, 30}, {20, 30}};
        vector<int> resultados_dist;
        for (const auto& par : pares_dist) {
            if (par.first > n || par.second > n) continue;
            resultados_dist.push_back(distancia(*grafo, par.first, par.second));
        }

        cout << "Analise concluida. Gerando relatorio...\n";

        ofstream out(arquivo_saida);
        out << "=== RELATORIO DE ANALISE DE GRAFO ===\n\n";
        out << "Arquivo de entrada: " << caminho_arquivo << "\n";
        out << "Tipo: " << (eh_direcionado ? "Direcionado" : "Nao-Direcionado") << "\n";
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
            if(bfs_pais_resultados.count(s)) {
                out << "  [BFS] Pais de (10, 20, 30): (" << bfs_pais_resultados[s][0] << ", " << bfs_pais_resultados[s][1] << ", " << bfs_pais_resultados[s][2] << ")\n";
            }
            if(dfs_pais_resultados.count(s)) {
                out << "  [DFS] Pais de (10, 20, 30): (" << dfs_pais_resultados[s][0] << ", " << dfs_pais_resultados[s][1] << ", " << dfs_pais_resultados[s][2] << ")\n";
            }
        }
        out << "\n";
        out << "--- Item 4.5: Distancia entre Pares ---\n";
        if(resultados_dist.size() >= 3) {
            out << "  Distancia (10, 20): " << resultados_dist[0] << "\n";
            out << "  Distancia (10, 30): " << resultados_dist[1] << "\n";
            out << "  Distancia (20, 30): " << resultados_dist[2] << "\n";
        }
        
        out.close();
        cout << "\nRelatorio salvo em '" << arquivo_saida << "'.\n";

    } 
    else if (modo == "dijkstra") {
        
        map<string, int> nome_para_id;
        vector<string> id_para_nome;
        string caminho_mapa = "";
        try {
            std::filesystem::path p_grafo(caminho_arquivo);
            std::filesystem::path p_mapa = p_grafo.parent_path() / "rede_colaboracao_vertices.txt";
            if (std::filesystem::exists(p_mapa)) {
                caminho_mapa = p_mapa.string();
            } else {
                p_mapa = "rede_colaboracao_vertices.txt";
                if (std::filesystem::exists(p_mapa)) {
                    caminho_mapa = p_mapa.string();
                }
            }
            if (!caminho_mapa.empty()) {
                cout << "Lendo mapeamento de vertices de '" << caminho_mapa << "'...\n";
                tie(nome_para_id, id_para_nome) = ler_mapeamento_vertices(caminho_mapa, n);
            } else {
                cout << "Aviso: Nao foi possivel encontrar 'rede_colaboracao_vertices.txt'.\n";
                cout << "       O Estudo de Caso 3.3 (Rede de Colaboracao) sera pulado.\n";
            }
        } catch (const std::filesystem::filesystem_error& e) {
            cerr << "Erro de filesystem ao tentar localizar 'rede_colaboracao_vertices.txt': " << e.what() << "\n";
            caminho_mapa = ""; 
        }
        
        if (tem_peso_negativo) {
            string msg_erro = "AVISO: O grafo possui arestas com peso negativo.\nA biblioteca ainda nao implementa caminhos minimos com pesos negativos.";
            cout << "\n------------------------------------------------------------------\n" << msg_erro << "\n------------------------------------------------------------------" << endl;
            cerr << msg_erro << endl; 
            return 0;
        }

        const string arquivo_saida_dijkstra = "Dijkstra_" + nome_base_arquivo + ".txt";
        ofstream out_dijkstra(arquivo_saida_dijkstra);
        if (!out_dijkstra) {
            cerr << "Erro: Nao foi possivel criar o arquivo de saida: " << arquivo_saida_dijkstra << endl;
            return 1;
        }

        cout << "Iniciando estudos de caso Dijkstra... (Resultados serao salvos em " << arquivo_saida_dijkstra << ")" << endl;
        
        out_dijkstra << "=== RELATORIO DE ANALISE DIJKSTRA ===\n\n";
        out_dijkstra << "Arquivo de entrada: " << caminho_arquivo << "\n";
        out_dijkstra << "Tipo: " << (eh_direcionado ? "Direcionado" : "Nao-Direcionado") << "\n";
        out_dijkstra << "Representacao utilizada: " << tipo_representacao << "\n";
        out_dijkstra << "Vertices: " << n << ", Arestas: " << m << "\n";

        out_dijkstra << "\n--- Estudo de Caso 3.1: Distancias e Caminhos (IDs) ---" << endl;
        int origem_fixa = 10;
        if (origem_fixa > n) {
             out_dijkstra << "Vertice de origem 10 nao existe no grafo." << endl;
        } else {
            out_dijkstra << "Origem: Vertice " << origem_fixa << endl;
            
            FilaComHeap fila_heap_estudo1; 
            auto [dist, pai] = dijkstra(*grafo, origem_fixa, fila_heap_estudo1);
            
            vector<int> destinos = {20, 30, 40, 50, 60};
            const double INFINITO = std::numeric_limits<double>::infinity();

            RelatorioUtils::imprimir_borda_tabela_dijkstra(out_dijkstra);
            RelatorioUtils::imprimir_linha_tabela(out_dijkstra, "Destino", "Distancia", "Caminho (por IDs)");
            RelatorioUtils::imprimir_borda_tabela_dijkstra(out_dijkstra);
            
            out_dijkstra << fixed << setprecision(4);
            for (int d : destinos) {
                if (d > n) continue;
                
                string dist_str = (dist[d] == INFINITO) ? "INF" : to_string(dist[d]);
                
                string caminho_str_completo = (dist[d] == INFINITO) ? "Inalcançavel" : RelatorioUtils::get_caminho_str(origem_fixa, d, pai);
                
                string caminho_str_tabela;
                bool foi_truncado = false;
                const int LIMITE_LARGURA_CAMINHO_IDS = 33; 

                if (caminho_str_completo.length() > LIMITE_LARGURA_CAMINHO_IDS) {
                    caminho_str_tabela = caminho_str_completo.substr(0, LIMITE_LARGURA_CAMINHO_IDS - 3) + "...";
                    foi_truncado = true;
                } else {
                    caminho_str_tabela = caminho_str_completo;
                }

                RelatorioUtils::imprimir_linha_tabela(out_dijkstra, to_string(d), dist_str, caminho_str_tabela);

                if (foi_truncado) {
                    out_dijkstra << "   |    ... Caminho Completo: " << caminho_str_completo << endl;
                }
            }
            RelatorioUtils::imprimir_borda_tabela_dijkstra(out_dijkstra);
        }
        

        int k_benchmark = 100; 
        string tipo_fila = "ambos"; 
        
        // Parsing inteligente dos argumentos opcionais do Dijkstra
        // Ignora "direcionado" pois já foi tratado
        if (argc >= 5) {
             // Tenta encontrar k_benchmark e tipo_fila nos argumentos restantes
             // que NÃO sejam "direcionado"
             vector<string> args_restantes;
             for(int i=4; i<argc; ++i) {
                 if(string(argv[i]) != "direcionado") {
                     args_restantes.push_back(argv[i]);
                 }
             }
             
             if(args_restantes.size() >= 1) {
                 try {
                    k_benchmark = stoi(args_restantes[0]);
                 } catch(...) {
                    tipo_fila = args_restantes[0];
                 }
             }
             if(args_restantes.size() >= 2) {
                 tipo_fila = args_restantes[1];
             }
        }

        if (tipo_fila != "vetor" && tipo_fila != "heap" && tipo_fila != "ambos") {
            cerr << "Aviso: Tipo de fila '" << tipo_fila << "' invalido. Usando default 'ambos'.\n";
            tipo_fila = "ambos";
        }

        out_dijkstra << "\n--- Estudo de Caso 3.2: Benchmark de Performance (Dijkstra) ---" << endl;
        out_dijkstra << "Executando " << k_benchmark << " vezes (Tipo Fila: " << tipo_fila << ")..." << endl;

        vector<int> origens_aleatorias;
        std::mt19937 gen(42); 
        std::uniform_int_distribution<> dis(1, n);
        for (int i = 0; i < k_benchmark; ++i) {
            origens_aleatorias.push_back(dis(gen));
        }

        long long t_total_vetor_ms = 0;
        double t_medio_vetor = 0.0;
        long long t_total_heap_ms = 0;
        double t_medio_heap = 0.0;

        if (tipo_fila == "vetor" || tipo_fila == "ambos") {
            cout << "Executando benchmark com Vetor...\n";
            auto start_vetor = chrono::high_resolution_clock::now();
            for (int s : origens_aleatorias) {
                FilaComVetor fila_vetor;
                volatile auto res = dijkstra(*grafo, s, fila_vetor); 
            }
            auto end_vetor = chrono::high_resolution_clock::now();
            t_total_vetor_ms = chrono::duration_cast<chrono::milliseconds>(end_vetor - start_vetor).count();
            t_medio_vetor = (k_benchmark > 0) ? (double)t_total_vetor_ms / k_benchmark : 0.0;
        }

        if (tipo_fila == "heap" || tipo_fila == "ambos") {
            cout << "Executando benchmark com Heap...\n";
            auto start_heap = chrono::high_resolution_clock::now();
            for (int s : origens_aleatorias) {
                FilaComHeap fila_heap;
                volatile auto res = dijkstra(*grafo, s, fila_heap);
            }
            auto end_heap = chrono::high_resolution_clock::now();
            t_total_heap_ms = chrono::duration_cast<chrono::milliseconds>(end_heap - start_heap).count();
            t_medio_heap = (k_benchmark > 0) ? (double)t_total_heap_ms / k_benchmark : 0.0;
        }

        RelatorioUtils::imprimir_tabela_benchmark(out_dijkstra, 
            t_total_vetor_ms, t_medio_vetor, 
            t_total_heap_ms, t_medio_heap);


        out_dijkstra << "\n--- Estudo de Caso 3.3: Rede de Colaboracao ---" << endl;
        
        if (caminho_mapa.empty() || nome_para_id.empty()) {
            out_dijkstra << "AVISO: Arquivo 'rede_colaboracao_vertices.txt' nao foi encontrado ou esta vazio.\n";
            out_dijkstra << "       Este estudo de caso foi pulado." << endl;
        
        } else {
            const string NOME_ORIGEM = "Edsger W. Dijkstra";
            const vector<string> NOMES_DESTINO = {
                "Alan M. Turing",     
                "Joseph B. Kruskal",  
                "Jon M. Kleinberg",   
                "Éva Tardos",         
                "Daniel Figueiredo"   
            };
            const double INFINITO = std::numeric_limits<double>::infinity();

            int id_origem = -1;
            if (nome_para_id.count(NOME_ORIGEM)) {
                id_origem = nome_para_id.at(NOME_ORIGEM);
            }

            if (id_origem == -1) {
                out_dijkstra << "Erro: Pesquisador de origem '" << NOME_ORIGEM << "' nao encontrado no mapeamento." << endl;
            
            } else {
                out_dijkstra << "Origem: " << NOME_ORIGEM << " (ID: " << id_origem << ")\n\n";
                
                FilaComHeap fila_heap_estudo3; 
                auto [dist, pai] = dijkstra(*grafo, id_origem, fila_heap_estudo3);

                RelatorioUtils::imprimir_borda_tabela_nomes(out_dijkstra);
                RelatorioUtils::imprimir_linha_tabela_nomes(out_dijkstra, "Destino", "Distancia", "Caminho Minimo");
                RelatorioUtils::imprimir_borda_tabela_nomes(out_dijkstra);
                
                out_dijkstra << fixed << setprecision(4);

                for (const string& nome_d : NOMES_DESTINO) {
                    int id_d = -1;
                    if (nome_para_id.count(nome_d)) {
                        id_d = nome_para_id.at(nome_d);
                    }

                    string dist_str;
                    string caminho_str_completo;
                    
                    if (id_d == -1) {
                        dist_str = "N/A";
                        caminho_str_completo = "Pesquisador '" + nome_d + "' nao encontrado no mapa.";
                    } else if (dist[id_d] == INFINITO) {
                        dist_str = "INF";
                        caminho_str_completo = "Inalcançavel";
                    } else {
                        dist_str = to_string(dist[id_d]);
                        caminho_str_completo = RelatorioUtils::get_caminho_str_com_nomes(id_origem, id_d, pai, id_para_nome);
                    }
                    
                    string caminho_str_tabela; 
                    bool foi_truncado = false;
                    const int LIMITE_LARGURA_CAMINHO = 48;

                    if (caminho_str_completo.length() > LIMITE_LARGURA_CAMINHO) {
                        caminho_str_tabela = caminho_str_completo.substr(0, LIMITE_LARGURA_CAMINHO - 3) + "...";
                        foi_truncado = true;
                    } else {
                        caminho_str_tabela = caminho_str_completo;
                    }
                    
                    RelatorioUtils::imprimir_linha_tabela_nomes(out_dijkstra, nome_d, dist_str, caminho_str_tabela);

                    if (foi_truncado) {
                        out_dijkstra << " |    ... Caminho Completo: " << caminho_str_completo << endl;
                    }
                }
                RelatorioUtils::imprimir_borda_tabela_nomes(out_dijkstra);
            }
        }
        
        out_dijkstra.close();
        cout << "\nRelatorio Dijkstra salvo com sucesso em '" << arquivo_saida_dijkstra << "'." << endl;

    }
    // ======================================================================
    // --- MODO 3: BELLMAN-FORD (Tarefa 3) ---
    // ======================================================================
    else if (modo == "bellman-ford") {
        
        // O enunciado pede:
        // 1. Distância de 10, 20, 30 para 100.
        // 2. Tempo médio de 10 rodadas.
        
        int destino_t = 100;
        if (n < 100) {
            cout << "Aviso: O grafo tem menos de 100 vertices. Mudando destino para " << n << ".\n";
            destino_t = n;
        }

        cout << "Executando Bellman-Ford (Todos -> " << destino_t << ")...\n";
        cout << "Rodando 10 vezes para calcular tempo medio...\n";

        // Variáveis para armazenar resultado
        vector<double> dist_final;
        vector<int> pai_final;
        bool ciclo_detectado = false;

        // Benchmark
        auto start = chrono::high_resolution_clock::now();
        for (int i = 0; i < 10; ++i) {
            auto resultado = bellman_ford(*grafo, destino_t);
            
            // Na última rodada, salvamos os dados para imprimir
            if (i == 9) {
                ciclo_detectado = resultado.first;
                dist_final = resultado.second.first;
                pai_final = resultado.second.second;
            }
        }
        auto end = chrono::high_resolution_clock::now();
        
        // Cálculos de tempo
        long long tempo_total_ms = chrono::duration_cast<chrono::milliseconds>(end - start).count();
        double tempo_medio_s = (double)tempo_total_ms / 10.0 / 1000.0; // ms -> s

        // --- SAÍDA DOS DADOS ---

        if (ciclo_detectado) {
            cout << "\nERRO CRITICO: Ciclo negativo detectado no grafo!\n";
            cout << "As distancias calculadas nao sao confiaveis.\n";
        } else {
            cout << "\n--- Estudo 3.1: Distancias para o vertice " << destino_t << " ---\n";
            RelatorioUtils::imprimir_borda_tabela_dijkstra(cout); // Reusando a borda
            RelatorioUtils::imprimir_linha_tabela(cout, "Origem", "Distancia", "Status");
            RelatorioUtils::imprimir_borda_tabela_dijkstra(cout);
            
            vector<int> origens = {10, 20, 30};
            const double INF = std::numeric_limits<double>::infinity();

            for (int u : origens) {
                if (u > n) continue;

                string d_str = (dist_final[u] == INF) ? "INF" : to_string(dist_final[u]);
                string status = (dist_final[u] == INF) ? "Inalcancavel" : "OK";
                
                RelatorioUtils::imprimir_linha_tabela(cout, to_string(u), d_str, status);
            }
            RelatorioUtils::imprimir_borda_tabela_dijkstra(cout);
        }

        cout << "\n--- Estudo 3.2: Performance ---\n";
        cout << "Tempo medio de execucao: " << fixed << setprecision(6) << tempo_medio_s << " segundos.\n";
        
        // Salvar num arquivo também, se quiser
        string arquivo_saida = "Relatorio_BellmanFord_" + nome_base_arquivo + ".txt";
        ofstream out_bf(arquivo_saida);
        out_bf << "Resultados Bellman-Ford (" << nome_base_arquivo << ")\n";
        out_bf << "Destino: " << destino_t << "\n";
        out_bf << "Tempo Medio: " << tempo_medio_s << " s\n";
        if(!ciclo_detectado) {
             out_bf << "Distancia 10->" << destino_t << ": " << dist_final[10] << "\n";
             out_bf << "Distancia 20->" << destino_t << ": " << dist_final[20] << "\n";
             out_bf << "Distancia 30->" << destino_t << ": " << dist_final[30] << "\n";
        }
        out_bf.close();
        cout << "Relatorio salvo em " << arquivo_saida << endl;
    }

    else {
        cerr << "Erro: Modo '" << modo << "' invalido. Use 'relatorio' ou 'dijkstra'.\n";
        return 1;
    }

    return 0;
}