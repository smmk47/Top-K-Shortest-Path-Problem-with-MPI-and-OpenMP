#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>

#define MAX_NODES 265214
#define K 3

typedef struct Edge {
    int vertex;
    int cost;
} Edge;

typedef struct AdjListNode {
    Edge edge;
    struct AdjListNode* next;
} AdjListNode;

typedef struct {
    AdjListNode* head;
} AdjList;

typedef struct Graph {
    int V;
    AdjList* array;
} Graph;

typedef struct HeapNode {
    int vertex;
    int dist;
} HeapNode;

typedef struct {
    HeapNode* nodes;
    int size;
    int capacity;
} MinHeap;

Graph* createGraph(int V) {
    Graph* graph = (Graph*) malloc(sizeof(Graph));
    graph->V = V;
    graph->array = (AdjList*) malloc(V * sizeof(AdjList));
    for (int i = 0; i < V; i++) {
        graph->array[i].head = NULL;
    }
    return graph;
}

void addEdge(Graph* graph, int src, int dest, int cost) {
    AdjListNode* newNode = (AdjListNode*) malloc(sizeof(AdjListNode));
    newNode->edge.vertex = dest;
    newNode->edge.cost = cost;
    newNode->next = graph->array[src].head;
    graph->array[src].head = newNode;
}

void freeGraph(Graph* graph) {
    for (int i = 0; i < graph->V; i++) {
        AdjListNode* node = graph->array[i].head;
        while (node) {
            AdjListNode* temp = node;
            node = node->next;
            free(temp);
        }
    }
    free(graph->array);
    free(graph);
}

MinHeap* createMinHeap(int capacity) {
    MinHeap* minHeap = (MinHeap*) malloc(sizeof(MinHeap));
    minHeap->nodes = (HeapNode*) malloc(capacity * sizeof(HeapNode));
    minHeap->size = 0;
    minHeap->capacity = capacity;
    return minHeap;
}

void insertMinHeap(MinHeap* minHeap, int v, int dist) {
    if (minHeap->size == minHeap->capacity) {
        return;
    }
    int i = minHeap->size++;
    minHeap->nodes[i].vertex = v;
    minHeap->nodes[i].dist = dist;
    while (i && minHeap->nodes[(i - 1) / 2].dist > minHeap->nodes[i].dist) {
        HeapNode tmp = minHeap->nodes[i];
        minHeap->nodes[i] = minHeap->nodes[(i - 1) / 2];
        minHeap->nodes[(i - 1) / 2] = tmp;
        i = (i - 1) / 2;
    }
}

HeapNode extractMin(MinHeap* minHeap) {
    if (minHeap->size <= 0) {
        return (HeapNode){-1, INT_MAX};
    }
    HeapNode root = minHeap->nodes[0];
    minHeap->nodes[0] = minHeap->nodes[--minHeap->size];
    int i = 0;
    while ((2 * i + 1) < minHeap->size) {
        int left = 2 * i + 1;
        int right = 2 * i + 2;
        int smallest = left;
        if (right < minHeap->size && minHeap->nodes[right].dist < minHeap->nodes[left].dist) {
            smallest = right;
        }
        if (minHeap->nodes[i].dist <= minHeap->nodes[smallest].dist) break;
        HeapNode tmp = minHeap->nodes[i];
        minHeap->nodes[i] = minHeap->nodes[smallest];
        minHeap->nodes[smallest] = tmp;
        i = smallest;
    }
    return root;
}

void readGraphFromFile(Graph* graph, const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Could not open file: %s\n", filename);
        return;
    }
    char line[256];
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '#') continue;
        int src, dest;
        if (sscanf(line, "%d\t%d", &src, &dest) == 2) {
            addEdge(graph, src, dest, 1);
        }
    }
    fclose(file);
}

char* findKShortestPaths(Graph* graph, int src, int dest) {
    int V = graph->V;
    int dis[V][K];
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < K; j++) {
            dis[i][j] = INT_MAX;
        }
    }

    MinHeap* minHeap = createMinHeap(V * K);
    insertMinHeap(minHeap, src, 0);
    dis[src][0] = 0;

    while (minHeap->size != 0) {
        HeapNode heapNode = extractMin(minHeap);
        int u = heapNode.vertex;

        for (AdjListNode* crawl = graph->array[u].head; crawl != NULL; crawl = crawl->next) {
            int v = crawl->edge.vertex;
            int weight = crawl->edge.cost;
            if (dis[v][K-1] > dis[u][0] + weight) {
                dis[v][K-1] = dis[u][0] + weight;
                for (int i = K-1; i > 0 && dis[v][i] < dis[v][i-1]; i--) {
                    int temp = dis[v][i];
                    dis[v][i] = dis[v][i-1];
                    dis[v][i-1] = temp;
                }
                insertMinHeap(minHeap, v, dis[v][K-1]);
            }
        }
    }

    char* paths = (char*)malloc(256 * K * sizeof(char));
    char temp[256];
    paths[0] = '\0'; 
    printf("Source: %d, Destination: %d\n", src, dest);
    for (int i = 0; i < K; i++) {
        if (dis[dest][i] == INT_MAX) {
            snprintf(temp, sizeof(temp), "Path %d: Infinity\n", i + 1);
        } else {
            snprintf(temp, sizeof(temp), "Path %d: %d\n", i + 1, dis[dest][i]);
        }
        strcat(paths, temp);
    }

    free(minHeap->nodes);
    free(minHeap);

    return paths;
}

int getMaxNodeID(Graph* graph) {
    int max = 0;
    for (int i = 0; i < graph->V; i++) {
        AdjListNode* node = graph->array[i].head;
        while (node) {
            if (node->edge.vertex > max) {
                max = node->edge.vertex;
            }
            node = node->next;
        }
    }
    return max;
}

void generateRandomPairs(int* sources, int* destinations, int num_pairs, int max_node) {
    srand(time(NULL));
    for (int i = 0; i < num_pairs; i++) {
        sources[i] = rand() % max_node;
        destinations[i] = rand() % max_node;
    }
}


int main() {
    clock_t start_time, end_time;
    double cpu_time_used;

    start_time = clock();

    // Load graph from file
    const char* filename = "Email-EuAll.txt";
    Graph* graph = createGraph(MAX_NODES);
    readGraphFromFile(graph, filename);

    int max_node = getMaxNodeID(graph);

    const int num_pairs = 10;
    int sources[num_pairs];
    int destinations[num_pairs];

    generateRandomPairs(sources, destinations, num_pairs, max_node);

    // 10 pairs
    for (int i = 0; i < num_pairs; i++) {
        char* paths = findKShortestPaths(graph, sources[i], destinations[i]);
        printf("%s", paths); 
        free(paths); 
    }

    freeGraph(graph);

    end_time = clock();
    cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

    printf("Execution time: %f seconds\n", cpu_time_used);

    return 0;
}