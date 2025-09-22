def ler_grafo(caminho):
    with open(caminho, 'r') as f:
        n = int(f.readline().strip())
        adj = [[] for _ in range(n + 1)]
        m = 0
        for linha in f:
            s = linha.strip()
            if not s:
                continue
            u, v = map(int, s.split())
            adj[u].append(v)
            adj[v].append(u)
            m += 1
    return n, m, adj

def matriz_de_adjacencia(n, adj):
    M = [[0]*n for _ in range(n)]
    for u in range(1, n+1):
        for v in adj[u]:
            M[u-1][v-1] = 1
            M[v-1][u-1] = 1
    return M

def insertion_sort(v):  # in-place
    for i in range(1, len(v)):
        x = v[i]
        j = i-1
        while j >= 0 and v[j] > x:
            v[j+1] = v[j]
            j -= 1
        v[j+1] = x

def graus_vertices(adj):
    return [len(adj[v]) for v in range(1, len(adj))]

def estatisticas_grau(graus):
    n = len(graus)
    if n == 0:
        return 0, 0, 0.0, 0
    gs = graus[:]
    insertion_sort(gs)
    mid = n // 2
    med = (gs[mid-1] + gs[mid]) / 2.0 if n % 2 == 0 else gs[mid]
    soma = 0.0
    for g in graus:
        soma += g
    return gs[0], gs[-1], soma/n, med

def bfs_lista(adj, s):
    n = len(adj) - 1
    pai = [-1]*(n+1)
    nivel = [-1]*(n+1)
    fila = [0]*(n+5); head = 0; tail = 0
    def push(x):
        nonlocal tail
        fila[tail] = x; tail += 1
    def pop():
        nonlocal head
        x = fila[head]; head += 1; return x
    nivel[s] = 0; push(s)
    while head < tail:
        u = pop()
        for v in adj[u]:
            if nivel[v] == -1:
                nivel[v] = nivel[u] + 1
                pai[v] = u
                push(v)
    return pai, nivel

def bfs_matriz(M, s):
    n = len(M)
    pai = [-1]*(n+1)
    nivel = [-1]*(n+1)
    fila = [0]*(n+5); head = 0; tail = 0
    def push(x):
        nonlocal tail
        fila[tail] = x; tail += 1
    def pop():
        nonlocal head
        x = fila[head]; head += 1; return x
    nivel[s] = 0; push(s)
    while head < tail:
        u = pop()
        row = M[u-1]
        i = 0
        while i < n:
            if row[i] == 1:
                v = i+1
                if nivel[v] == -1:
                    nivel[v] = nivel[u] + 1
                    pai[v] = u
                    push(v)
            i += 1
    return pai, nivel

def dfs_lista(adj, s):
    n = len(adj) - 1
    pai = [-1]*(n+1)
    nivel = [-1]*(n+1)
    cor = [0]*(n+1)  # 0 branco, 1 cinza, 2 preto
    pilha = [(s, 0)]
    nivel[s] = 0; cor[s] = 1
    while pilha:
        u, it = pilha.pop()
        if it < len(adj[u]):
            v = adj[u][it]
            pilha.append((u, it+1))
            if cor[v] == 0:
                cor[v] = 1
                pai[v] = u
                nivel[v] = nivel[u] + 1
                pilha.append((v, 0))
        else:
            cor[u] = 2
    return pai, nivel

def dfs_matriz(M, s):
    n = len(M)
    pai = [-1]*(n+1)
    nivel = [-1]*(n+1)
    cor = [0]*(n+1)
    pilha = [(s, 0)]
    nivel[s] = 0; cor[s] = 1
    while pilha:
        u, i = pilha.pop()
        # avança até achar vizinho branco
        while i < n and (M[u-1][i] == 0 or cor[i+1] != 0):
            i += 1
        if i < n:
            v = i+1
            pilha.append((u, i+1))
            cor[v] = 1
            pai[v] = u
            nivel[v] = nivel[u] + 1
            pilha.append((v, 0))
        else:
            cor[u] = 2
    return pai, nivel

def componentes(adj):
    n = len(adj) - 1
    visto = [0]*(n+1)
    comps = []
    for s in range(1, n+1):
        if visto[s] == 0:
            _, niv = bfs_lista(adj, s)
            comp = []
            for v in range(1, n+1):
                if niv[v] != -1 and visto[v] == 0:
                    visto[v] = 1
                    comp.append(v)
            comps.append(comp)
    # ordenar por tamanho desc (insertion sort)
    for i in range(1, len(comps)):
        x = comps[i]; j = i-1
        while j >= 0 and len(comps[j]) < len(x):
            comps[j+1] = comps[j]; j -= 1
        comps[j+1] = x
    return comps

def distancia(adj, u, v):
    _, nivel = bfs_lista(adj, u)
    return nivel[v]

def diametro_exato(adj):
    n = len(adj) - 1
    diam = -1
    for s in range(1, n+1):
        _, nivel = bfs_lista(adj, s)
        mx = -1
        for v in range(1, n+1):
            if nivel[v] > mx:
                mx = nivel[v]
        if mx > diam:
            diam = mx
    return diam

def farthest_from(adj, s):
    _, nivel = bfs_lista(adj, s)
    far = s
    for v in range(1, len(adj)):
        if nivel[v] > nivel[far]:
            far = v
    return far, nivel[far]

def diametro_aproximado(adj):
    # 2-BFS: pega arbitrário a -> b mais distante; depois de b -> c mais distante; usa dist(b,c)
    if len(adj) <= 2: return 0
    a = 1
    b, _ = farthest_from(adj, a)
    c, dist = farthest_from(adj, b)
    return dist

def pausa(msg):
    print(msg)
    _ = input()  # só pra pausar sem importações

def estudos_de_caso(caminho='grafo_1.txt', saida='ResultadoGrafo.txt'):
    n, m, adj = ler_grafo(caminho)
    M = matriz_de_adjacencia(n, adj)

    # Item 1 (memória): pausas para você medir no SO
    pausa("REPRESENTAÇÃO: LISTA carregada. Meça a MEMÓRIA do processo e pressione ENTER...")
    pausa("REPRESENTAÇÃO: MATRIZ carregada. Meça a MEMÓRIA do processo e pressione ENTER...")

    # Itens 2 e 3 (tempos): 100 BFS/DFS em cada representação
    # Você mede o tempo externamente (cronômetro/`time`) entre as mensagens PRONTO/FIM.
    seeds = []
    # sementes 1..n repetindo até 100 execuções
    if n > 0:
        v = 1
        while len(seeds) < 100:
            seeds.append(v)
            v += 1
            if v > n: v = 1

    print("Item 2: 100 BFS na LISTA - PRONTO. Aperte ENTER para iniciar, cronometre externamente.")
    _ = input()
    for s in seeds:
        bfs_lista(adj, s)
    print("Item 2: 100 BFS na LISTA - FIM.")

    print("Item 2: 100 BFS na MATRIZ - PRONTO. Aperte ENTER para iniciar.")
    _ = input()
    for s in seeds:
        bfs_matriz(M, s)
    print("Item 2: 100 BFS na MATRIZ - FIM.")

    print("Item 3: 100 DFS na LISTA - PRONTO. Aperte ENTER para iniciar.")
    _ = input()
    for s in seeds:
        dfs_lista(adj, s)
    print("Item 3: 100 DFS na LISTA - FIM.")

    print("Item 3: 100 DFS na MATRIZ - PRONTO. Aperte ENTER para iniciar.")
    _ = input()
    for s in seeds:
        dfs_matriz(M, s)
    print("Item 3: 100 DFS na MATRIZ - FIM.")

    # Item 4: pais de 10,20,30 iniciando em 1,2,3
    consultar = [10, 20, 30]
    fontes = [1, 2, 3]
    pais_bfs = {}; pais_dfs = {}
    for s in fontes:
        pb, _ = bfs_lista(adj, s)
        pd, _ = dfs_lista(adj, s)
        pais_bfs[s] = [pb[v] if v < len(pb) else None for v in consultar]
        pais_dfs[s] = [pd[v] if v < len(pd) else None for v in consultar]

    # Item 5: distâncias
    pares = [(10,20),(10,30),(20,30)]
    dists = []
    for (a,b) in pares:
        if a <= n and b <= n:
            dists.append(distancia(adj, a, b))
        else:
            dists.append(None)

    # Item 6: componentes
    comps = componentes(adj)
    qtd = len(comps)
    tamanhos = [len(c) for c in comps]
    maior = tamanhos[0] if tamanhos else 0
    menor = tamanhos[-1] if tamanhos else 0

    # Item 7: diâmetro (exato) e aproximado
    diam = diametro_exato(adj)
    diam_aprox = diametro_aproximado(adj)

    # Estatísticas de grau (pra relatório completo)
    graus = graus_vertices(adj)
    gmin, gmax, gmed, gmediana = estatisticas_grau(graus)

    with open(saida, 'w') as out:
        out.write("=== ESTUDOS DE CASO ===\n")
        out.write(f"Vértices: {n}\nArestas: {m}\n")
        out.write(f"Grau mín/máx/médio/mediana: {gmin} / {gmax} / {gmed:.6f} / {gmediana}\n\n")

        out.write("Item 1 (Memória): meça manualmente nas pausas do programa.\n\n")

        out.write("Item 2 (BFS) & Item 3 (DFS): cronometre externamente entre PRONTO e FIM.\n\n")

        out.write("Item 4 - Pais de {10,20,30} nas árvores geradas:\n")
        for s in fontes:
            out.write(f"  BFS a partir de {s}: pais(10,20,30) = {pais_bfs[s]}\n")
        for s in fontes:
            out.write(f"  DFS a partir de {s}: pais(10,20,30) = {pais_dfs[s]}\n")
        out.write("\n")

        out.write("Item 5 - Distâncias (10,20), (10,30), (20,30):\n")
        out.write(f"  {pares[0]} = {dists[0]}\n")
        out.write(f"  {pares[1]} = {dists[1]}\n")
        out.write(f"  {pares[2]} = {dists[2]}\n\n")

        out.write("Item 6 - Componentes Conexas:\n")
        out.write(f"  Quantidade: {qtd}\n")
        out.write(f"  Maior: {maior}  |  Menor: {menor}\n")
        out.write(f"  Tamanhos (desc): {tamanhos}\n\n")

        out.write("Item 7 - Diâmetro:\n")
        out.write(f"  Exato: {diam}\n")
        out.write(f"  Aproximado (2-BFS): {diam_aprox}\n")

    print("Relatório salvo em", saida)


# Execute:
# Coloque seu arquivo 'grafo_1.txt' na mesma pasta e rode o script.
estudos_de_caso()
