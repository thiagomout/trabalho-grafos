from os import write


def grafo():
    arestas = 0
    grau = 0
    try:
        with open('grafo_1.txt', 'r') as grafo1:
            vertices = grafo1.readline()
            vertices = int(vertices) ## Pega a linha 1 do arquivo e declara como vértice.
            grau = [0]*(vertices+1) ## Grau de cada vértice 0

            for linha in grafo1: ## Lendo arestas de u pra v, map fazendo com que a leitura seja feita em todo o txt
                u,v = map(int, linha.split())
                arestas += 1
                grau[u] += 1
                grau[v] += 1

    except FileNotFoundError:
        print("Não tem grafo.")
    ####################
    graus_cresc = sorted(grau)
    mid = vertices // 2
    if vertices % 2 == 0:            ##caso seja par
        mediana = (graus_cresc[mid-1]+graus_cresc[mid])/2
    else:
        mediana = graus_cresc[mid]    ##caso seja ímpar
    print("----------")
    print("Vertices:",vertices)
    grau = grau[1::]
    print("Graus",grau)
    gmax = (max(grau))
    gmin = (min(grau))
    gemd = (sum(grau)/vertices)
    narestas = len(grau)

    result = "ResultadoGrafo.txt"
    try:
        with open(result,'w') as novoGrafo:
            novoGrafo.write("Número de Vertices: "+ str(vertices)+"\n")
            novoGrafo.write("Grau máximo: "+ str(gmax) +"\n")
            novoGrafo.write("Grau mínimo: " + str(gmin) + "\n")
            novoGrafo.write("Mediana: "+ str(mediana)+"\n")
            novoGrafo.write("Grau Médio: "+ str(gemd) + "\n")
            novoGrafo.write("Numero de Arestas: "+ str(narestas)+"\n")

    except:
        print("Não há caminho pro arquivo")

grafo()