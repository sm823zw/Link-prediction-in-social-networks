#include <stdio.h>
#include <stdlib.h>
#include "17990_linked_list.h"

// Data structure to store vertices connecting an edge and the corresponding score
struct Edge{
    int first;
    int second;
    int steps;
    float score;
};

// Function to create an Edge
struct Edge* newEdge(int first, int second, float score){
    struct Edge* edge = malloc(sizeof(struct Edge));
    edge->first = first;
    edge->second = second;
    edge->steps = 0;
    edge->score = score;
    return edge;
}

// Function to create a specialized Edge
struct Edge* newEdge2(int first, int second, int steps, float score){
    struct Edge* edge = malloc(sizeof(struct Edge));
    edge->first = first;
    edge->second = second;
    edge->steps = steps;
    edge->score = score;
    return edge;
}

// Comparator function to compare two edges while sorting
int comparator_function(const void *v1, const void *v2){
    struct Edge* e1 = *(struct Edge**) v1;
    struct Edge* e2 = *(struct Edge**) v2;
    if((float)e1->score < (float)e2->score){
        return +1;
    }
    else if((float)e1->score > (float)e2->score){
        return -1;
    }
    else{
        if(e1->first > e2->first){
            return +1;
        }
        else if(e1->first < e2->first){
            return -1;
        }
        else{
            if(e1->second >= e2->second){
                return +1;
            }
            else{
                return -1;
            }
        }
    }
    return 0;
}

// Data structure to store adjacency list
struct AdjList{
    int data;
    int index;
    struct Node* neighbours;
    struct AdjList* nextVertex;
};

// Create a new adjacency list for a vertex
struct AdjList* newAdjList(int data, int idx){
    struct AdjList* adjlist = malloc(sizeof(struct AdjList));
    adjlist->data = data;
    adjlist->index = idx;
    adjlist->neighbours = NULL;
    adjlist->nextVertex = NULL;
    return adjlist;
}

// Function to fetch the index of element needed to be accessed in a 3D matrix
int idx(int i, int j, int k, int N_VERTICES){
  return i + (j * (N_VERTICES)) + (k * (N_VERTICES) * (N_VERTICES));
}

// Function that returns the list of neighbours for a given vertex
struct Node* neighbour_list(struct AdjList* adj, int x){
    struct AdjList* curr = adj;
    while(curr){
        if(curr->data == x){
            return curr->neighbours;
        }
        curr = curr->nextVertex;
    }
}

// Function that returns the index assigned to a vertex
int vertex_index(struct AdjList* adj, int x){
    struct AdjList* curr = adj;
    while(curr){
        if(curr->data == x){
            return curr->index;
        }
        curr = curr->nextVertex;
    }
}

// Function that checks whether two vertices are neighbours or not
int neighbour_check(struct AdjList* adj, int v1, int v2){
    struct AdjList* x = adj;
    while(x){
        if(x->data == v1){
            return member(x->neighbours, v2);
            break;
        }
        x = x->nextVertex;
    }
    return 0;
}

// Function that adds an edge to the adjacency list
void addEdge(int src, int des, struct AdjList* adj){
    struct AdjList* curr = adj;
    while(curr){
        if(curr->data == src){
            struct Node* x = curr->neighbours;
            if(!x){
                curr->neighbours = newNode(des);
                return;
            }
            while(x->next){
                x = x->next;
            }
            x->next = newNode(des);
            return;
        }
        curr = curr->nextVertex;
    }
}

// Function that prints the adjacency list
void print_adj_list(struct AdjList* adj){
    FILE *adjfile = fopen("adjList.txt", "w");
    struct AdjList* curr = adj;
    while(curr){
        fprintf(adjfile, "%d -> ", curr->data);
        struct Node* x = curr->neighbours;
        while(x->next){
            fprintf(adjfile, "%d -> ", x->data);
            x = x->next;
        }
        fprintf(adjfile, "%d\n", x->data);
        curr = curr->nextVertex;
    }
    fclose(adjfile);
}

// Function to compute Jaccard score for all non-existent edges and sort them
void jaccard(struct Node* vertex_list, struct AdjList* adj, int N_VERTICES, int N_EDGES, int K){
    int total_edges = (N_VERTICES*(N_VERTICES-1)/2) - N_EDGES;
    struct Edge* edges[total_edges];
    int cnt = 0;
    struct Node* x = vertex_list;
    while(x){
        struct Node* y = x->next;
        while(y){
            if(!neighbour_check(adj, x->data, y->data)){
                struct Node* x_neighbour = neighbour_list(adj, x->data);
                struct Node* y_neighbour = neighbour_list(adj, y->data);
                int a = set_union(x_neighbour, y_neighbour);
                int b = set_intersection(x_neighbour, y_neighbour);
                float score = (float)b/(float)a;
                edges[cnt] = newEdge(x->data, y->data, score);
                cnt++;
            }
            y = y->next;
        }
        x = x->next;
    }    
    qsort(edges, total_edges, sizeof(struct Edge*), comparator_function);
    FILE *write_fptr = fopen("Jaccard.txt", "w");
    for(int i=0;i<K;i++){
        fprintf(write_fptr, "%d %d %f\n", edges[i]->first, edges[i]->second, edges[i]->score);
        free(edges[i]);
    }
    fclose(write_fptr);
}

// Function to get count of all paths of lengths upto 6
void multiply_matrix(struct AdjList* adj, int N_VERTICES, long A[][N_VERTICES], long* path_count){
    long a_2[N_VERTICES][N_VERTICES];
    for(int i=0;i<N_VERTICES;i++){
        for(int j=0;j<N_VERTICES;j++){
            a_2[i][j] = 0;
            path_count[idx(i, j, 2, N_VERTICES)] = 0;
            for(int k=0;k<N_VERTICES;k++){
                a_2[i][j] += A[i][k]*A[k][j];
            }
            path_count[idx(i, j, 2, N_VERTICES)] = a_2[i][j];
        }
    }
    long a_3[N_VERTICES][N_VERTICES];
    for(int i=0;i<N_VERTICES;i++){
        for(int j=0;j<N_VERTICES;j++){
            a_3[i][j] = 0;
            path_count[idx(i, j, 3, N_VERTICES)] = 0;
            for(int k=0;k<N_VERTICES;k++){
                a_3[i][j] += a_2[i][k]*A[k][j];
            }
            path_count[idx(i, j, 3, N_VERTICES)] = a_3[i][j];
        }
    }
    for(int i=0;i<N_VERTICES;i++){
        for(int j=0;j<N_VERTICES;j++){
            path_count[idx(i, j, 4, N_VERTICES)] = 0;
            path_count[idx(i, j, 5, N_VERTICES)] = 0;
            path_count[idx(i, j, 6, N_VERTICES)] = 0;
            for(int k=0;k<N_VERTICES;k++){
                path_count[idx(i, j, 4, N_VERTICES)] += a_2[i][k]*a_2[k][j];
                path_count[idx(i, j, 5, N_VERTICES)] += a_3[i][k]*a_2[k][j];
                path_count[idx(i, j, 6, N_VERTICES)] += a_3[i][k]*a_3[k][j];
            }
        }
    }
}

// Function to compute the Katz score of an edge
float katz_score(int src, int des, float beta, int N_VERTICES, struct AdjList* adj, long* path_count){
    float score = 0.0;
    float beta_p = beta*beta;
    int src_idx = vertex_index(adj, src);
    int des_idx = vertex_index(adj, des);
    for(int len=2;len<7;len++){
        score += beta_p*path_count[idx(src_idx, des_idx, len, N_VERTICES)];
        beta_p = beta_p*beta;
    }
    return score;
}

// Function to compute Katz score for all non-existent edges and sort them
void katz(int N_VERTICES, int N_EDGES, struct Node* vertex_list, struct AdjList* adj, long A[][N_VERTICES], long* path_count, float beta, int K){
    multiply_matrix(adj, N_VERTICES, A, path_count);
    int total_edges = (N_VERTICES*(N_VERTICES-1)/2) - N_EDGES;
    struct Edge* edges[total_edges];
    int cnt = 0;
    struct Node* x = vertex_list;
    while(x){
        struct Node* y = x->next;
        while(y){
            if(!neighbour_check(adj, x->data, y->data)){
                edges[cnt] = newEdge(x->data, y->data, katz_score(x->data, y->data, beta, N_VERTICES, adj, path_count));
                cnt++;
            }
            y = y->next;
        }
        x = x->next;
    }
    qsort(edges, total_edges, sizeof(struct Edge*), comparator_function);
    FILE *write_fptr = fopen("Katz.txt", "w");
    for(int i=0;i<K;i++){
        fprintf(write_fptr, "%d %d %f\n", edges[i]->first, edges[i]->second, edges[i]->score);
        free(edges[i]);
    }
    fclose(write_fptr);
}

// Function to get probability to go from one vertex to other in given steps
float get_prob(int src, int des, int m, int N_VERTICES, float* prob_tran){
    float p = 0.0;
    for(int i=0;i<N_VERTICES;i++){
        if(i != des){
            p += prob_tran[idx(src, i, 1, N_VERTICES)] * prob_tran[idx(i, des, m-1, N_VERTICES)];
        }
    }
    return p;
}

// Function to compute probability transition matrices till given number of steps
void compute_prob_tran_matrix(int steps, int N_VERTICES, struct Node* vertex_list, float* prob_tran, struct AdjList* adj){
    for(int i=0;i<N_VERTICES;i++){
        for(int j=0;j<N_VERTICES;j++){
            for(int k=0;k<steps+1;k++){
                prob_tran[idx(i, j, k, N_VERTICES)] = 0.0;
            }
        }
    }
    struct Node* x = vertex_list;
    while(x){
        struct Node* nl = neighbour_list(adj, x->data);
        int neighbour_count = length_ll(nl);
        struct Node* y = nl;
        while(y){
            prob_tran[idx(vertex_index(adj, x->data), vertex_index(adj, y->data), 1, N_VERTICES)] = (float)1/(float)neighbour_count;
            y = y->next;
        }
        x = x->next;
    }
    for(int m=2;m<steps+1;m++){
        for(int i=0;i<N_VERTICES;i++){
            for(int j=0;j<N_VERTICES;j++){
                float p = get_prob(i, j, m, N_VERTICES, prob_tran);
                prob_tran[idx(i, j, m, N_VERTICES)] = p;
            }
        }
    }
}

// Function to calculate the expected time of random walk
float calc_expectation(int src, int des, int steps, int N_VERTICES, float* prob_tran){
    float expectation = 0.0;
    for(int i=2;i<steps+1;i++){
        expectation += i * prob_tran[idx(src, des, i, N_VERTICES)];
    }
    return expectation;
}

// Function to compute Commute time for all non-existent edges and sort them
void commuteTime(int steps, int N_VERTICES, int N_EDGES, struct Node* vertex_list, float* prob_tran, struct AdjList* adj, int K){
    compute_prob_tran_matrix(steps, N_VERTICES, vertex_list, prob_tran, adj);
    int total_edges = (N_VERTICES*(N_VERTICES-1)/2) - N_EDGES;
    struct Edge* edges[total_edges];
    int cnt = 0;
    struct Node* x = vertex_list;
    while(x){
        struct Node* y = x->next;
        while(y){
            if(!neighbour_check(adj, x->data, y->data)){
                int x_idx = vertex_index(adj, x->data);
                int y_idx = vertex_index(adj, y->data);
                float exp = calc_expectation(x_idx, y_idx, steps, N_VERTICES, prob_tran);
                exp += calc_expectation(y_idx, x_idx, steps, N_VERTICES, prob_tran);
                edges[cnt] = newEdge(x->data, y->data, -exp);
                cnt++;
            }
            y = y->next;
        }
        x = x->next;
    }
    qsort(edges, total_edges, sizeof(struct Edge*), comparator_function);
    FILE *write_fptr = fopen("HittingTime.txt", "w");
    for(int i=0;i<K;i++){
        fprintf(write_fptr, "%d %d %f\n", edges[i]->first, edges[i]->second, edges[i]->score);
        free(edges[i]);
    }
    fclose(write_fptr);
}

// Function to compute accurate Commute time for all the non-existent edges, compute the steps taken and sort them
void commuteTimeAccurate(int steps, int N_VERTICES, int N_EDGES, struct Node* vertex_list, float* prob_tran, struct AdjList* adj, int K){
    compute_prob_tran_matrix(steps, N_VERTICES, vertex_list, prob_tran, adj);
    int total_edges = (N_VERTICES*(N_VERTICES-1)/2) - N_EDGES;
    struct Edge* edges[total_edges];
    int cnt = 0;
    struct Node* x = vertex_list;
    while(x){
        struct Node* y = x->next;
        while(y){
            if(!neighbour_check(adj, x->data, y->data)){
                int x_idx = vertex_index(adj, x->data);
                int y_idx = vertex_index(adj, y->data);
                float exp = 0.0;
                float prev = 0.0;
                int steps_accurate = 0;
                for(int i=2;i<steps+1;i++){
                    float abc = i * prob_tran[idx(x_idx, y_idx, i, N_VERTICES)];
                    abc += i * prob_tran[idx(y_idx, x_idx, i, N_VERTICES)];
                    if(prev>0.0 && abc >= 0.01){
                        exp = calc_expectation(x_idx, y_idx, i, N_VERTICES, prob_tran);
                        exp += calc_expectation(y_idx, x_idx, i, N_VERTICES, prob_tran);
                        steps_accurate = i;
                        break;
                    }
                    prev = abc;
                }
                edges[cnt] = newEdge2(x->data, y->data, steps_accurate, -exp);
                cnt++;
            }
            y = y->next;
        }
        x = x->next;
    }
    qsort(edges, total_edges, sizeof(struct Edge*), comparator_function);
    FILE *write_fptr = fopen("HittingTimeAccurate.txt", "w");
    for(int i=0;i<K;i++){
        fprintf(write_fptr, "%d %d %f %d\n", edges[i]->first, edges[i]->second, edges[i]->score, edges[i]->steps);
        free(edges[i]);
    }
    fclose(write_fptr);
}

// Function that returns the rooted page rank given an edge
float calc_pageRank(int v1, int v2, int N_VERTICES, float* prob_tran, float alpha, int itr){
    float transition_mat[N_VERTICES][N_VERTICES];
    for(int i=0;i<N_VERTICES;i++){
        for(int j=0;j<N_VERTICES;j++){
            transition_mat[i][j] = (float)prob_tran[idx(i, j, 1, N_VERTICES)];
        }
    }
    for(int i=0;i<N_VERTICES;i++){
        transition_mat[v2][i] *= (1-alpha);
        transition_mat[i][v1] = alpha;
    }
    float eig_vec[N_VERTICES];
    for(int i=0;i<N_VERTICES;i++){
        eig_vec[i] = 1.0;
    }
    while(itr>0){
        float xyz[N_VERTICES];
        float max = 0;
        for(int i=0;i<N_VERTICES;i++){
            xyz[i] = 0;
            for(int j=0;j<N_VERTICES;j++){
                xyz[i] += eig_vec[j]*transition_mat[i][j];
            }
            if(max < xyz[i]){
                max = xyz[i];
            }
        }
        for(int i=0;i<N_VERTICES;i++){
            eig_vec[i] = (float)xyz[i]/(float)max;            
        }
        itr--;
    }
    float sum = 0;
    for(int i=0;i<N_VERTICES;i++){
        sum += eig_vec[i];
    }
    for(int i=0;i<N_VERTICES;i++){
        eig_vec[i] /= sum;
    }
    return eig_vec[v2];
}

// Function to compute rooted Page Rank for all the non-existent edges and sort them
void pageRank(int N_VERTICES, int N_EDGES, struct Node* vertex_list, float* prob_tran, struct AdjList* adj, float alpha, int itr, int K){
    int total_edges = (N_VERTICES*(N_VERTICES-1)/2) - N_EDGES;
    struct Edge* edges[total_edges];
    int cnt = 0;
    struct Node* x = vertex_list;
    while(x){
        struct Node* y = x->next;
        while(y){
            if(!neighbour_check(adj, x->data, y->data)){
                int x_idx = vertex_index(adj, x->data);
                int y_idx = vertex_index(adj, y->data);
                float score = calc_pageRank(x_idx, y_idx, N_VERTICES, prob_tran, alpha, itr);
                score += calc_pageRank(y_idx, x_idx, N_VERTICES, prob_tran, alpha, itr);
                score /= 2;
                edges[cnt] = newEdge(x->data, y->data, score);
                cnt++;
            }
            y = y->next;
        }
        x = x->next;
    }
    qsort(edges, total_edges, sizeof(struct Edge*), comparator_function);
    FILE *write_fptr = fopen("PageRank.txt", "w");
    for(int i=0;i<K;i++){
        fprintf(write_fptr, "%d %d %f\n", edges[i]->first, edges[i]->second, edges[i]->score);
        free(edges[i]);
    }
    fclose(write_fptr);
}

// Function that returns adjacency list of the graph
struct AdjList* compute_adjacency_list(struct Node* vertex_list, int N_EDGES, FILE* fptr){
    if(fptr == NULL){
        printf("Error opening file");
        exit(1);
    }
    struct AdjList* adj = newAdjList(-1, -1);
    struct AdjList* curr = adj;
    struct Node* y = vertex_list;
    int idx = 0;
    while(y){
        curr->nextVertex = newAdjList(y->data, idx);
        curr = curr->nextVertex;
        y = y->next;
        idx++;
    }
    adj = adj->nextVertex;
    int v1, v2, w1;
    int c = 0;
    while(c < N_EDGES){
        fscanf(fptr, "%d %d %d", &v1, &v2, &w1);
        // fscanf(fptr, "%d %d", &v1, &v2);
        addEdge(v1, v2, adj);
        addEdge(v2, v1, adj);
        c++;
    }
    return adj;
}

// Function that returns adjacency matrix of the graph (required and used only for Katz score computation)
void compute_adjacency_matrix(struct AdjList* adj, int N_EDGES, int N_VERTICES, long A[][N_VERTICES], FILE *fptr){
    if(fptr == NULL){
        printf("Error opening file");
        exit(1);
    }
    for(int i=0;i<N_VERTICES;i++){
        for(int j=0;j<N_VERTICES;j++){
            A[i][j] = 0;
        }
    }
    int v1, v2, w1;
    int c = 0;
    while(c < N_EDGES){
        fscanf(fptr, "%d %d %d", &v1, &v2, &w1);
        // fscanf(fptr, "%d %d", &v1, &v2);
        A[vertex_index(adj, v1)][vertex_index(adj, v2)] = 1;
        A[vertex_index(adj, v2)][vertex_index(adj, v1)] = 1;
        c++;
    }
}

// Function to count number of edges in a graph
int count_edges(FILE* fptr){
    if(fptr == NULL){
        printf("Error opening file");
        exit(1);
    }
    int N_EDGES = 0;
    for(char i=getc(fptr); i!=EOF;i=getc(fptr)){
        if(i == '\n'){
            N_EDGES++;
        }
    }
    fclose(fptr);
    return N_EDGES;
}

// Function that returns a list of vertices
struct Node* list_of_vertex(int N_EDGES, FILE* fptr){
    if(fptr == NULL){
        printf("Error opening file");
        exit(1);
    }
    int v1, v2, w1;
    int c = 0;
    struct Node* vertex_list = newNode(-1);
    while(c < N_EDGES){
        fscanf(fptr, "%d %d %d", &v1, &v2, &w1);
        // fscanf(fptr, "%d %d", &v1, &v2);
        vertex_list = add_vertex(vertex_list, v1);
        vertex_list = add_vertex(vertex_list, v2);
        c++;
    }
    return vertex_list->next;
}

int main(){

    char file_name[] = "../contact-high-school-proj-graph.txt";

    FILE* fptr = fopen(file_name, "r");
    int N_EDGES = count_edges(fptr);
    fclose(fptr);

    fptr = fopen(file_name, "r");
    struct Node* vertex_list = list_of_vertex(N_EDGES, fptr);
    int N_VERTICES = length_ll(vertex_list);
    fclose(fptr);

    fptr = fopen(file_name, "r");
    struct AdjList* adj = compute_adjacency_list(vertex_list, N_EDGES, fptr);
    fclose(fptr);
    print_adj_list(adj);
    printf("Adjacency list printed. Check adjList.txt.\n");

    int K = 2;
    printf("Enter Value for K -> ");
    scanf("%d", &K);

    // Compute Jaccard Score
    jaccard(vertex_list, adj, N_VERTICES, N_EDGES, K);
    printf("Jaccard scores for top %d edges printed. Check Jaccard.txt.\n", K);

    // Compute Katz Score
    fptr = fopen(file_name, "r");
    long A[N_VERTICES][N_VERTICES];
    compute_adjacency_matrix(adj, N_EDGES, N_VERTICES, A, fptr);
    fclose(fptr);
    long* path_count = malloc(sizeof(path_count) * (N_VERTICES) * (N_VERTICES) * 7);
    float beta = 0.1;
    katz(N_VERTICES, N_EDGES, vertex_list, adj, A, path_count, beta, K);
    free(path_count);
    printf("Katz scores for top %d edges printed. Check Katz.txt.\n", K);

    // Initialization for the next two parts
    int steps = 50;
    float* prob_tran = malloc(sizeof(prob_tran) * (N_VERTICES) * (N_VERTICES) * (steps+1));
    
    // Compute Commute Time
    steps = 6;
    commuteTime(steps, N_VERTICES, N_EDGES, vertex_list, prob_tran, adj, K);
    printf("Commute time for top %d edges printed. Check HittingTime.txt.\n", K);

    // Compute Commute Time Accurately 
    steps = 50; 
    commuteTimeAccurate(steps, N_VERTICES, N_EDGES, vertex_list, prob_tran, adj, K);
    printf("Commute time accurately for top %d edges printed. Check HittingTimeAccurate.txt.\n", K);

    // Rooted Page Rank
    float alpha = 0.2;
    int itr = 2; // Number of iterations of the power method. Takes a lot of time for this graph. Can be increased for smaller graphs 
    pageRank(N_VERTICES, N_EDGES, vertex_list, prob_tran, adj, alpha, itr, K);
    printf("Rooted Page Rank for top %d edges printed. Check PageRank.txt.\n", K);
    free(prob_tran);


    return 0;
}
