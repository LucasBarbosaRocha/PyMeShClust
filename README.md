# Clustering MeanShift (PyMeShClust)
* Nome: Lucas Barbosa Rocha
* Disciplina: Inteligência Artificial
* Trabalho: Implementar um clustering para sequências de DNA utilizando MeanShift.
* Contato: lucas.lb.rocha@gmail.com

Projeto irmão: **PyKmClust** (mesmo problema, resolvido com K-means). Os dois foram feitos para serem comparados.

## A ideia
Baseado no artigo *MeShClust: an intelligent tool for clustering DNA sequences*. Os autores implementaram em C++; esta versão é em Python com numpy.

Agrupar sequências de DNA alinhando todas contra todas é caro: são O(n²) alinhamentos. A alternativa é comparar sequências **sem alinhar** (*alignment-free*), pela composição de k-mers:

1. **Sequência → vetor.** Cada sequência vira um vetor com a frequência de cada k-mer (palavra de tamanho k). Para k=3 são 64 posições: AAA, AAC, ..., TTT. O vetor é normalizado (soma 1).
   ```
   ACGTAC  ->  ACG, CGT, GTA, TAC  ->  [0, ..., 0.25 (ACG), ..., 0.25 (CGT), ...]
   ```
2. **Similaridade.** A similaridade entre duas sequências é a interseção dos vetores, `soma(min(a, b))`, entre 0 (nada em comum) e 1 (mesma composição).
3. **Fase 1, mean shift guloso.** Começa com uma sequência como centro:
   * todas as sequências livres com similaridade ≥ limiar entram no cluster;
   * **mean shift:** a média dos vetores dos membros é o *centro sintético*; o novo centro é o membro real mais parecido com essa média. O centro "anda" em direção à região mais densa;
   * repete a busca com o novo centro até ninguém mais entrar;
   * o próximo cluster começa pela sequência livre mais parecida com o último centro.
4. **Fase 2, junção.** Clusters cujos centros têm similaridade ≥ limiar são juntados, do maior para o menor, e o centro é recalculado. Repete até estabilizar.
5. **Saída.** Um arquivo `.clstr` no formato do CD-HIT; o centro de cada cluster é marcado com `*`.

A diferença central para o PyKmClust: aqui **você não define o número de clusters**, e sim o quão parecidas as sequências de um cluster precisam ser. O número de clusters é consequência.

### Diferença para o MeShClust do artigo
O MeShClust treina um GLM que prevê a **identidade de alinhamento** a partir de estatísticas de k-mers e aplica o limiar sobre essa previsão. Aqui o limiar é aplicado direto na similaridade de k-mers, então `-s 0.95` **não** equivale a 95% de identidade no CD-HIT. Implementar o GLM exigiria alinhar pares de sequências para gerar os dados de treino.

## Como usar

### Instalação
A única dependência é o `numpy`. Em Ubuntu/Debian o `pip` não instala no Python do sistema (erro `externally-managed-environment`), então use uma das opções:

```bash
# Opção 1: numpy do sistema
sudo apt install python3-numpy

# Opção 2: ambiente virtual (precisa do pacote python3-venv)
sudo apt install python3.12-venv
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
source .venv/bin/activate
```

### Rodar o clustering
```bash
python3 main.py sequencias.fasta -s 0.95
python3 main.py hbv.fasta -s 0.8 --kmer 6 -o saida.clstr
```

| Parâmetro | Padrão | Descrição |
|---|---|---|
| `entrada` | `sequencias.fasta` | arquivo FASTA (sequência em uma ou várias linhas) |
| `-s` | 0.95 | limiar de similaridade entre 0 e 1 |
| `--kmer` | 3 | tamanho do k-mer |
| `--sem-juncao` | | pula a fase 2 |
| `-o` | `output<s>.clstr` | arquivo de saída |

A similaridade cai quando o k-mer cresce (k-mers longos se repetem menos entre sequências). **Com k-mer maior, use limiar menor.**

### Avaliar um resultado
`avaliar.py` compara um `.clstr` com um rótulo verdadeiro extraído do cabeçalho das sequências:

```bash
python3 avaliar.py output95.clstr --rotulo especie
```

| Rótulo | Exemplo de cabeçalho | Classe |
|---|---|---|
| `genotipo` | `HBV genotype C DNA, complete genome` | C |
| `regiao` | `HBV genotype F S gene ..., partial cds` | gene S parcial |
| `especie` | `Klebsiella pneumoniae strain SIKP041` | Klebsiella pneumoniae |
| `genero` | `Klebsiella pneumoniae strain SIKP041` | Klebsiella |

Sequências sem rótulo (os reads `ERR599000...` do `sequencias.fasta` e 4 sequências do HBV sem genótipo) são agrupadas normalmente, mas ficam fora da avaliação.

Métricas (1 = perfeito):
* **ARI** (Adjusted Rand Index): concordância entre os pares de sequências, corrigida pelo acaso. ~0 = agrupamento aleatório. É a métrica principal.
* **NMI**: informação mútua normalizada entre clusters e classes.
* **Pureza**: fração das sequências que estão no cluster dominado pela sua classe. Sobe sozinha quando há muitos clusters pequenos (com um cluster por sequência a pureza é 1), então olhe sempre junto com o ARI.

ARI e NMI foram conferidos contra o scikit-learn.

### Reproduzir os experimentos
```bash
python3 experimentos.py   # gera RESULTADOS.md e resultados.csv
```

## Resultados
Tabelas completas em [RESULTADOS.md](RESULTADOS.md) e [resultados.csv](resultados.csv), com limiares de 0.5 a 0.95. Valores em ARI.

**Dados:**
* `hbv.fasta`: 16 sequências de HBV, com 7 genótipos e 4 regiões do genoma (genoma completo, gene S, gene P, P/preS/S).
* `sequencias.fasta`: 500 sequências, sendo 237 contigs de 13 espécies/8 gêneros de bactérias e 263 reads curtos (~100nt) sem rótulo.

### Efeito do tamanho do k-mer (melhor limiar para cada k-mer)
| k-mer | HBV região | HBV genótipo | gênero, todas | gênero, só rotuladas | matriz p/ 1 milhão de seqs |
|---|---|---|---|---|---|
| 3 | **0.966** (s=0.95) | 0.065 (s=0.9) | **0.197** (s=0.8) | **0.194** (s=0.8) | 256 MB |
| 4 | **0.966** (s=0.9) | 0.316 (s=0.95) | 0.128 (s=0.7) | 0.129 (s=0.7) | 1 GB |
| 5 | **0.966** (s=0.8) | 0.316 (s=0.9) | 0.090 (s=0.5) | 0.090 (s=0.5) | 4 GB |
| 6 | 0.801 (s=0.5) | **0.552** (s=0.8) | 0.054 (s=0.5) | 0.054 (s=0.5) | 16 GB |

* **O melhor limiar cai quando o k-mer cresce**, como esperado: 0.95 no k-mer 3 e 0.8 no k-mer 5 para a mesma região do HBV.
* **Nas bactérias, com k-mer 5 e 6, mesmo s=0.5 fragmenta demais** (95 e 149 clusters para 237 sequências). A grade de limiares deveria descer abaixo de 0.5 para esses casos.
* **O genótipo do HBV só aparece com k-mer 6.** K-mers longos capturam diferenças mais finas de sequência.
* **Os reads curtos quase não afetam o PyMeShClust** (as colunas "todas" e "só rotuladas" são praticamente iguais). Reads que não se parecem com nada formam clusters próprios em vez de puxar os outros, ao contrário do k-means.

### Comparação PyKmClust × PyMeShClust
Melhor configuração de cada método (entre parênteses: k-mer e parâmetro).

| Tarefa | PyKmClust (k-means) | PyMeShClust |
|---|---|---|
| HBV região | 0.801 (3, k=4) | **0.966** (3, s=0.95, 6 clusters) |
| HBV genótipo | 0.300 (4, k=7) | **0.552** (6, s=0.8, 13 clusters) |
| Bactérias, espécie (todas) | **0.223** (4, k=8) | 0.198 (3, s=0.8, 247 clusters) |
| Bactérias, gênero (só rotuladas) | **0.329** (5, k=8) | 0.194 (3, s=0.8, 25 clusters) |

## Conclusões
1. **A representação importa mais que o algoritmo.** A diferença entre k-means e MeShClust é menor do que a diferença causada pelo tamanho do k-mer ou pela presença de reads curtos. Escolher bem como a sequência vira vetor vem antes de escolher o clustering.
2. **Nenhum dos dois ganha sempre, porque cada um supõe uma forma diferente para os grupos.**
   * O **k-means** força exatamente k grupos e divide o espaço em regiões. Vai melhor quando a classe é espalhada (contigs de bactérias vindos de regiões diferentes do genoma). Precisa saber o k e é sensível a outliers.
   * O **MeShClust** usa um limiar igual para todos os clusters. Encontra bem grupos compactos (regiões do HBV), mas quebra classes espalhadas em muitos clusters pequenos. O ruído forma clusters próprios em vez de contaminar os outros.
   * A escolha depende do formato dos dados e de saber ou não quantos grupos existem.
3. **A composição de k-mers mede "que tipo de sequência é", não "quão parecida ela é".** Os dois métodos separam bem a região do genoma e mal o genótipo, que difere em ~8% dos nucleotídeos: mutações pontuais quase não mudam a frequência de k-mers curtos. É por isso que o MeShClust do artigo usa um GLM para converter estatísticas de k-mers em identidade de alinhamento.
4. **Os dados limitam as conclusões.** O HBV tem só 16 sequências. Nas bactérias, contigs de ~800nt de genes diferentes da mesma espécie têm composições diferentes, então o ARI baixo mistura falha do método com uma tarefa que esse sinal não resolve. Mais da metade do `sequencias.fasta` não tem rótulo.
5. **A avaliação com rótulo é indispensável.** A quantidade de clusters ou a pureza sozinhas enganam: o MeShClust com limiar 0.95 tem pureza 0.99 nas bactérias, mas ARI 0.05, porque os clusters são quase unitários.
6. **Uso prático.** Clustering por k-mers sem alinhamento é rápido e escala bem (200 mil sequências em segundos, com pouca memória). Serve para agrupamento grosso e pré-filtragem: separar tipos de sequência, reduzir redundância, fazer um primeiro corte antes de alinhar. Não substitui métodos baseados em identidade para separar variantes próximas, como genótipos ou cepas.

### Próximos passos
* Um dataset rotulado maior: centenas de genomas completos de HBV com genótipo conhecido (por exemplo, do HBVdb).
* Comparar com o CD-HIT nos mesmos dados, para medir quanto se perde sem alinhamento.
* Implementar o GLM do MeShClust (alinhando uma amostra de pares para treino) ou usar distâncias MinHash no estilo Mash.

## Arquivos
* `main.py`: linha de comando do clustering.
* `kmer.py`: leitura de FASTA, vetores de k-mers e escrita do `.clstr`.
* `meshclust.py`: as duas fases do algoritmo.
* `avaliar.py`: rótulos, ARI, NMI e pureza.
* `experimentos.py`: gera `RESULTADOS.md` e `resultados.csv`.
* `original/`: primeira versão do código (2020).

## Memória
Todas as sequências ficam em uma matriz numpy `(N, 4^k)` float32 e um vetor de rótulos. As similaridades são calculadas em blocos de 65536 linhas. Cada k-mer a mais multiplica a memória por 4:

| Sequências | k=3 | k=4 | k=5 | k=6 |
|---|---|---|---|---|
| 100 mil | 25 MB | 100 MB | 400 MB | 1.6 GB |
| 1 milhão | 256 MB | 1 GB | 4 GB | 16 GB |

Com 200 mil sequências de 800nt (arquivo de 167MB), k-mer 3 e limiar 0.9: ~130MB de pico de RSS e ~4s. O tempo cresce com N × quantidade de clusters, então muitos clusters pequenos deixam a execução mais lenta.
