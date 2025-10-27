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
            vizinhos_ids.reserve(adj[u].size()); // linha de otimização de memória
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


// -------------------- Leitura do grafo --------------------
tuple<int, long long, vector<vector<pair<int, double>>>, bool> ler_dados_grafo(const string& caminho) {
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
        adj[u].push_back({v, w});
        adj[v].push_back({u, w});
        m++;
    }
    return {n, m, adj, tem_peso_negativo};
}

// -------------------- Algoritmos Genéricos --------------------

class FilaDePrioridade {
public:
    virtual ~FilaDePrioridade() = default;

    virtual void inicializar(int num_vertices, int s) = 0;
    
    virtual bool esta_vazia() const = 0;
    
    // Retorna o ID do vértice com a menor distância e o remove da fila.
    virtual int extrair_min() = 0;
    
    // Atualiza a distância de um vértice na fila (se for menor).
    virtual void inserir_ou_atualizar(int vertice, double nova_distancia) = 0;

    virtual double get_distancia(int v) const = 0;
};

class FilaComVetor : public FilaDePrioridade {
private:
    vector<double> dist;
    vector<bool> visitado;
    int n;
    int restantes; // Contador de quantos vértices ainda faltam visitar

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

        // Loop O(V): Encontra o vértice não-visitado com menor distância
        for (int v = 1; v <= n; ++v) {
            if (!visitado[v] && dist[v] < min_dist) {
                min_dist = dist[v];
                min_v = v;
            }
        }

        if (min_v != -1) {
            visitado[min_v] = true; // "Remove" da fila
            restantes--;
        }
        return min_v; // Retorna -1 se todos os restantes forem infinitos
    }

    void inserir_ou_atualizar(int vertice, double nova_distancia) override {
        // alterar o valor no vetor de distâncias.
        dist[vertice] = nova_distancia;
    }
    double get_distancia(int v) const override {
        return dist[v];
    }
};

class FilaComHeap : public FilaDePrioridade {
private:
// O 'greater' o torna um min-heap (ordena pelo menor primeiro).
    using Par = pair<double, int>; // {distancia, vertice}
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
            // Verifica se é a distância mais atual
            if (d == dist[v]) {
                return v;
            }
        }
        return -1; // Fila vazia
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
        if (u == -1 || fila.get_distancia(u) == INFINITO) break; // Todos os restantes são infinitos

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

// -------------------- Utilitários de Saída/Relatório --------------------
namespace RelatorioUtils {

    /**
     * @brief Retorna uma string com o caminho (ex: "1 -> 5 -> 10").
     */
    string get_caminho_str(int origem, int destino, const vector<int>& pai) {
        // Se o destino não tem pai e não é a própria origem...
        if (pai[destino] == -1 && destino != origem) {
            return (destino == origem) ? to_string(origem) : "Nenhum";
        }

        stringstream ss;
        stack<int> caminho;
        int atual = destino;
        while (atual != -1) { // -1 é o 'pai' da origem
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

    /**
     * @brief Imprime uma linha de tabela formatada em um stream (cout ou file).
     * @param out O stream de saída (ex: cout ou um ofstream).
     */
    void imprimir_linha_tabela(ostream& out, const string& col1, const string& col2, const string& col3) {
        out << "  | " << left << setw(10) << col1
            << " | " << left << setw(18) << col2
            << " | " << left << setw(35) << col3 << " |" << endl;
    }

    /**
     * @brief Imprime o cabeçalho/rodapé da tabela de distâncias do Dijkstra.
     */
    void imprimir_borda_tabela_dijkstra(ostream& out) {
        out << "  +------------+--------------------+-------------------------------------+" << endl;
    }

    /**
     * @brief Imprime a tabela de resultados do benchmark de performance.
     */
    void imprimir_tabela_benchmark(ostream& out, 
                                   long long t_total_vetor_ms, double t_medio_vetor, 
                                   long long t_total_heap_ms, double t_medio_heap) 
    {
        out << "  +------------------+-------------------+------------------+" << endl;
        out << "  | Implementacao    | Tempo Total (ms)  | Tempo Medio (ms) |" << endl;
        out << "  +------------------+-------------------+------------------+" << endl;
        out << fixed << setprecision(4); // Garante a formatação correta dos 'doubles'
        out << "  | Vetor (O(V^2))   | " << left << setw(17) << t_total_vetor_ms 
             << " | " << left << setw(16) << t_medio_vetor << " |" << endl;
        out << "  | Heap (O(E log V))| " << left << setw(17) << t_total_heap_ms 
             << " | " << left << setw(16) << t_medio_heap << " |" << endl;
        out << "  +------------------+-------------------+------------------+" << endl;
    }

} // fim do namespace RelatorioUtils

// -------------------- Driver --------------------
int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // --- 1. Validação dos Argumentos ---
    if (argc < 4) {
        cerr << "Erro: Uso incorreto.\n";
        cerr << "Modo Relatorio: ./estudos <lista|matriz> <arquivo.txt> relatorio\n";
        cerr << "Modo Dijkstra:  ./estudos <lista|matriz> <arquivo.txt> dijkstra [k_benchmark]\n";
        return 1;
    }
    string tipo_representacao = argv[1];
    if (tipo_representacao != "lista" && tipo_representacao != "matriz") {
        cerr << "Erro: Representacao '" << tipo_representacao << "' invalida. Use 'lista' ou 'matriz'.\n";
        return 1;
    }
    string caminho_arquivo = argv[2];
    string modo = argv[3];

    // --- Extração do Nome Base do Arquivo ---
    string nome_base_arquivo;

    size_t pos_barra = caminho_arquivo.find_last_of("/\\");
    if (pos_barra != string::npos) {
        nome_base_arquivo = caminho_arquivo.substr(pos_barra + 1);
    } else {
        nome_base_arquivo = caminho_arquivo;
    }
    // Remove a extensão (ex: .txt)
    size_t pos_ponto = nome_base_arquivo.find_last_of(".");
    if (pos_ponto != string::npos) {
        nome_base_arquivo = nome_base_arquivo.substr(0, pos_ponto);
    }


    // --- 2. Leitura do Grafo ---
    cout << "Lendo dados do grafo de '" << caminho_arquivo << "'...\n";
    auto [n, m, adj_data, tem_peso_negativo] = ler_dados_grafo(caminho_arquivo);
    
    unique_ptr<Grafo> grafo;
    if (tipo_representacao == "lista") {
        grafo = make_unique<GrafoListaAdj>(n, m, move(adj_data));
        cout << "Grafo carregado com representacao em LISTA de adjacencia.\n";
    } else {
        grafo = make_unique<GrafoMatrizAdj>(n, m, adj_data);
        cout << "Grafo carregado com representacao em MATRIZ de adjacencia.\n";
    }
    cout << "Vertices: " << n << ", Arestas: " << m << endl;


    // ======================================================================
    // --- MODO 1: Trabalho1 (BFS/DFS, Diâmetro, Componentes) ---
    // ======================================================================
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
        // ... (toda a sua lógica de escrita de relatório original vai aqui) ...
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
        
        // ... (sua lógica de pais e distâncias) ...
        out << "--- Item 4.4: Pais dos Vertices 10, 20, 30 ---\n";
        for (int s : origens_pais) {
            out << "Buscas a partir do vertice " << s << ":\n";
            // Adicionando verificação de existência da chave
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
    // ======================================================================
    // --- MODO 2: ESTUDOS DE CASO DIJKSTRA ---
    // ======================================================================
    else if (modo == "dijkstra") {
        
        // --- 3. Verificação de Pesos Negativos ---
        if (tem_peso_negativo) {
            string msg_erro = "AVISO: O grafo possui arestas com peso negativo.\nA biblioteca ainda nao implementa caminhos minimos com pesos negativos.";
            cout << "\n------------------------------------------------------------------\n" << msg_erro << "\n------------------------------------------------------------------" << endl;
            cerr << msg_erro << endl; // Erros também vão para stderr
            return 0; // Encerra conforme pedido
        }

        // --- [MODIFICADO] ABRIR ARQUIVO DE SAÍDA PARA O DIJKSTRA ---
        const string arquivo_saida_dijkstra = "Dijkstra_" + nome_base_arquivo + ".txt";
        ofstream out_dijkstra(arquivo_saida_dijkstra);
        if (!out_dijkstra) {
            cerr << "Erro: Nao foi possivel criar o arquivo de saida: " << arquivo_saida_dijkstra << endl;
            return 1;
        }

        cout << "Iniciando estudos de caso Dijkstra... (Resultados serao salvos em " << arquivo_saida_dijkstra << ")" << endl;
        
        // Cabeçalho do arquivo
        out_dijkstra << "=== RELATORIO DE ANALISE DIJKSTRA ===\n\n";
        out_dijkstra << "Arquivo de entrada: " << caminho_arquivo << "\n";
        out_dijkstra << "Representacao utilizada: " << tipo_representacao << "\n";
        out_dijkstra << "Vertices: " << n << ", Arestas: " << m << "\n";

        // --- 4. Estudo de Caso 3.1: Distância (10 -> X) ---
        out_dijkstra << "\n--- Estudo de Caso 3.1: Distancias e Caminhos ---" << endl;
        int origem_fixa = 10;
        if (origem_fixa > n) {
             out_dijkstra << "Vertice de origem 10 nao existe no grafo." << endl;
        } else {
            out_dijkstra << "Origem: Vertice " << origem_fixa << endl;
            
            FilaComHeap fila_heap_estudo1; 
            auto [dist, pai] = dijkstra(*grafo, origem_fixa, fila_heap_estudo1);
            
            vector<int> destinos = {20, 30, 40, 50, 60};
            const double INFINITO = std::numeric_limits<double>::infinity();

            // Usando as funções do namespace RelatorioUtils
            RelatorioUtils::imprimir_borda_tabela_dijkstra(out_dijkstra);
            RelatorioUtils::imprimir_linha_tabela(out_dijkstra, "Destino", "Distancia", "Caminho");
            RelatorioUtils::imprimir_borda_tabela_dijkstra(out_dijkstra);
            
            out_dijkstra << fixed << setprecision(4); // Formatação para os pesos 'double'
            for (int d : destinos) {
                if (d > n) continue; // Ignora se o vértice não existe no grafo
                
                string dist_str = (dist[d] == INFINITO) ? "INF" : to_string(dist[d]);
                string caminho_str = (dist[d] == INFINITO) ? "Inalcançavel" : RelatorioUtils::get_caminho_str(origem_fixa, d, pai);
                RelatorioUtils::imprimir_linha_tabela(out_dijkstra, to_string(d), dist_str, caminho_str);
            }
            RelatorioUtils::imprimir_borda_tabela_dijkstra(out_dijkstra);
        }

        // --- 5. Estudo de Caso 3.2: Benchmark Dijkstra (Vetor vs Heap) ---
        int k_benchmark = 100; // Default K
        if (argc > 4) { // Se o usuário forneceu um K
            try {
                k_benchmark = stoi(argv[4]);
            } catch (...) {
                cerr << "Aviso: k_benchmark invalido. Usando default (100).\n";
                k_benchmark = 100;
            }
        }

        out_dijkstra << "\n--- Estudo de Caso 3.2: Benchmark de Performance (Dijkstra) ---" << endl;
        out_dijkstra << "Executando " << k_benchmark << " vezes com origens aleatorias..." << endl;

        // Gerar k vértices de origem aleatórios
        vector<int> origens_aleatorias;
        std::mt19937 gen(42); // Seed fixa para reprodutibilidade
        std::uniform_int_distribution<> dis(1, n);
        for (int i = 0; i < k_benchmark; ++i) {
            origens_aleatorias.push_back(dis(gen));
        }

        // Benchmark com Vetor
        auto start_vetor = chrono::high_resolution_clock::now();
        for (int s : origens_aleatorias) {
            FilaComVetor fila_vetor;
            volatile auto res = dijkstra(*grafo, s, fila_vetor); 
        }
        auto end_vetor = chrono::high_resolution_clock::now();
        auto t_total_vetor_ms = chrono::duration_cast<chrono::milliseconds>(end_vetor - start_vetor).count();
        double t_medio_vetor = (k_benchmark > 0) ? (double)t_total_vetor_ms / k_benchmark : 0.0;

        // Benchmark com Heap
        auto start_heap = chrono::high_resolution_clock::now();
        for (int s : origens_aleatorias) {
            FilaComHeap fila_heap;
            volatile auto res = dijkstra(*grafo, s, fila_heap);
        }
        auto end_heap = chrono::high_resolution_clock::now();
        auto t_total_heap_ms = chrono::duration_cast<chrono::milliseconds>(end_heap - start_heap).count();
        double t_medio_heap = (k_benchmark > 0) ? (double)t_total_heap_ms / k_benchmark : 0.0;

        // Tabela de Resultados do Benchmark
        RelatorioUtils::imprimir_tabela_benchmark(out_dijkstra, 
            t_total_vetor_ms, t_medio_vetor, 
            t_total_heap_ms, t_medio_heap);


        // --- 6. Estudo de Caso 3.3: Rede de Colaboração ---
        out_dijkstra << "\n--- Estudo de Caso 3.3: Rede de Colaboracao ---" << endl;
        out_dijkstra << "Este estudo de caso requer um grafo com identificadores (strings) para os vertices." << endl;
        out_dijkstra << "A biblioteca atual so suporta vertices indexados por numeros inteiros (1 a N)." << endl;
        out_dijkstra << "Para implementa-lo, seria necessario modificar a leitura do grafo para " << endl;
        out_dijkstra << "mapear nomes (ex: 'Edsger W. Dijkstra') para IDs (ex: 1) e vice-versa." << endl;
        
        out_dijkstra.close();
        cout << "\nRelatorio Dijkstra salvo com sucesso em '" << arquivo_saida_dijkstra << "'." << endl;

    }
    // ======================================================================
    // --- MODO INVÁLIDO ---
    // ======================================================================
    else {
        cerr << "Erro: Modo '" << modo << "' invalido. Use 'relatorio' ou 'dijkstra'.\n";
        return 1;
    }

    return 0;
}