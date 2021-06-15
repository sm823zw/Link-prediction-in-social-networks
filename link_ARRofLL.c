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

// Function to fetch the index of element needed to be accessed in a 3D matrix
int idx(int i, int j, int k, int N_VERTICES){
  return i + (j * (N_VERTICES+1)) + (k * (N_VERTICES+1) * (N_VERTICES+1));
}

// Function to add edges for Undirected graphs
void addEdgeUDG(int src, int dest, struct Node* adj[]){
    struct Node* curr = adj[src];
    while(curr->next){
        curr = curr->next;
    }
    curr->next = newNode(dest);
    curr = adj[dest];
    while(curr->next){
        curr = curr->next;
    }
    curr->next = newNode(src);
}

// Function to print adjacency list of a graph
void printGraph(int N_VERTICES, struct Node* adj[]){
    FILE *adjfile = fopen("adjList.txt", "w");

    for(int i=0;i<N_VERTICES+1;i++){
        struct Node* curr = adj[i];
        while(curr->next){
            fprintf(adjfile, "%d -> ", curr->data);
            curr = curr->next;
        }
        fprintf(adjfile, "%d\n", curr->data);
    }
    fclose(adjfile);
}

// Function to compute Jaccard score for all non-existent edges and sort them
void jaccard(int N_EDGES, int N_VERTICES, struct Node* adj[], int K){
    int total_edges = (N_VERTICES*(N_VERTICES-1)/2) - N_EDGES;
    struct Edge* edges[total_edges];
    int cnt = 0;
    for(int i=1;i<N_VERTICES + 1;i++){
        for(int j=i+1;j<N_VERTICES + 1;j++){
            int flag = 0;
            struct Node* curr = adj[i];
            while(curr){
                if(curr->data == j){
                    flag = 1;
                    break;
                }
                curr = curr->next;
            }
            if(flag == 0){
                int x = set_union(adj[i]->next, adj[j]->next);
                int y = set_intersection(adj[i]->next, adj[j]->next);
                float score = (float)y/(float)x;
                edges[cnt] = newEdge(i, j, score);
                cnt++;
            }
        }
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
void multiply_matrix(int N_VERTICES, long A[][N_VERTICES+1], long* path_count){
    long a_2[N_VERTICES+1][N_VERTICES+1];
    for(int i=0;i<N_VERTICES+1;i++){
        for(int j=0;j<N_VERTICES+1;j++){
            a_2[i][j] = 0;
            path_count[idx(i, j, 2, N_VERTICES)] = 0;
            for(int k=0;k<N_VERTICES+1;k++){
                a_2[i][j] += A[i][k]*A[k][j];
            }
            path_count[idx(i, j, 2, N_VERTICES)] = a_2[i][j];
        }
    }
    long a_3[N_VERTICES+1][N_VERTICES+1];
    for(int i=0;i<N_VERTICES+1;i++){
        for(int j=0;j<N_VERTICES+1;j++){
            a_3[i][j] = 0;
            path_count[idx(i, j, 3, N_VERTICES)] = 0;
            for(int k=0;k<N_VERTICES+1;k++){
                a_3[i][j] += a_2[i][k]*A[k][j];
            }
            path_count[idx(i, j, 3, N_VERTICES)] = a_3[i][j];
        }
    }
    for(int i=0;i<N_VERTICES+1;i++){
        for(int j=0;j<N_VERTICES+1;j++){
            path_count[idx(i, j, 4, N_VERTICES)] = 0;
            path_count[idx(i, j, 5, N_VERTICES)] = 0;
            path_count[idx(i, j, 6, N_VERTICES)] = 0;
            for(int k=0;k<N_VERTICES+1;k++){
                path_count[idx(i, j, 4, N_VERTICES)] += a_2[i][k]*a_2[k][j];
                path_count[idx(i, j, 5, N_VERTICES)] += a_3[i][k]*a_2[k][j];
                path_count[idx(i, j, 6, N_VERTICES)] += a_3[i][k]*a_3[k][j];
            }
        }
    }
}

// Function to compute the Katz score of an edge
float katz_score(int src, int des, float beta, int N_VERTICES, long* path_count){
    float score = 0.0;
    float beta_p = beta*beta;
    for(int len=2;len<7;len++){
        score += beta_p*path_count[idx(src, des, len, N_VERTICES)];
        beta_p = beta_p*beta;
    }
    return score;
}

// Function to compute Katz score for all non-existent edges and sort them
void katz(int N_VERTICES, int N_EDGES, struct Node* adj[], long A[][N_VERTICES+1], long* path_count, float beta, int K){
    multiply_matrix(N_VERTICES, A, path_count);
    int total_edges = (N_VERTICES*(N_VERTICES-1)/2) - N_EDGES;
    struct Edge* edges[total_edges];
    int cnt = 0;
    for(int i=1;i<N_VERTICES+1;i++){
        for(int j=i+1;j<N_VERTICES+1;j++){
            int flag = 0;
            struct Node* curr = adj[i];
            while(curr){
                if(curr->data == j){
                    flag = 1;
                    break;
                }
                curr = curr->next;
            }
            if(flag == 0){
                edges[cnt] = newEdge(i, j, katz_score(i, j, beta, N_VERTICES, path_count));
                cnt++;
            }
        }
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
    for(int i=1;i<N_VERTICES+1;i++){
        if(i != des){
            p += prob_tran[idx(src, i, 1, N_VERTICES)] * prob_tran[idx(i, des, m-1, N_VERTICES)];
        }
    }
    return p;
}

// Function to compute probability transition matrices till given number of steps
void compute_prob_tran_matrix(int steps, int N_VERTICES, float* prob_tran, struct Node* adj[]){
    for(int i=0;i<N_VERTICES+1;i++){
        for(int j=0;j<N_VERTICES+1;j++){
            for(int k=0;k<steps+1;k++){
                prob_tran[idx(i, j, k, N_VERTICES)] = 0.0;
            }
        }
    }
    for(int i=1;i<N_VERTICES+1;i++){
        int neigh_count = count_neighbours(adj[i]);
        struct Node* curr = adj[i]->next;
        while(curr){
            prob_tran[idx(i, curr->data, 1, N_VERTICES)] = (float)1/(float)neigh_count;
            curr = curr->next;
        }
    }
    for(int m=2;m<steps+1;m++){
        for(int i=1;i<N_VERTICES+1;i++){
            for(int j=1;j<N_VERTICES+1;j++){
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
void commuteTime(int steps, int N_VERTICES, int N_EDGES, float* prob_tran, struct Node* adj[], int K){
    compute_prob_tran_matrix(steps, N_VERTICES, prob_tran, adj);
    int total_edges = (N_VERTICES*(N_VERTICES-1)/2) - N_EDGES;
    struct Edge* edges[total_edges];
    int cnt = 0;
    for(int i=1;i<N_VERTICES + 1;i++){
        for(int j=i+1;j<N_VERTICES + 1;j++){
            int flag = 0;
            struct Node* curr = adj[i];
            while(curr){
                if(curr->data == j){
                    flag = 1;
                    break;
                }
                curr = curr->next;
            }
            if(flag == 0){
                float exp = calc_expectation(i, j, steps, N_VERTICES, prob_tran);
                exp += calc_expectation(j, i, steps, N_VERTICES, prob_tran);
                edges[cnt] = newEdge(i, j, -exp);
                cnt++;
            }
        }
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
void commuteTimeAccurate(int steps, int N_VERTICES, int N_EDGES, float* prob_tran, struct Node* adj[], int K){
    compute_prob_tran_matrix(steps, N_VERTICES, prob_tran, adj);
    int total_edges = (N_VERTICES*(N_VERTICES-1)/2) - N_EDGES;
    struct Edge* edges[total_edges];
    int cnt = 0;
    for(int i=1;i<N_VERTICES + 1;i++){
        for(int j=i+1;j<N_VERTICES + 1;j++){
            int flag = 0;
            struct Node* curr = adj[i];
            while(curr){
                if(curr->data == j){
                    flag = 1;
                    break;
                }
                curr = curr->next;
            }
            if(flag == 0){
                float exp;
                float prev = 0.0;
                int steps_accurate = 0;
                for(int k=2;k<steps+1;k++){
                    float abc = k*prob_tran[idx(i,j,k,N_VERTICES)] + k*prob_tran[idx(j,i,k,N_VERTICES)];
                    if(prev > 0.0 && abc >= 0.01){
                        exp = calc_expectation(i, j, k, N_VERTICES, prob_tran);
                        exp += calc_expectation(j, i, k, N_VERTICES, prob_tran);
                        steps_accurate = k;
                        break;
                    }
                    prev = abc;
                }
                edges[cnt] = newEdge2(i, j, steps_accurate, -exp);
                cnt++;
            }
        }
    }
    qsort(edges, total_edges, sizeof(struct Edge*), comparator_function);
    FILE *write_fptr = fopen("HittingTimeAccurate.txt", "w");
    for(int i=0;i<K;i++){
        fprintf(write_fptr, "%d %d %f %d\n", edges[i]->first, edges[i]->second, edges[i]->score, edges[i]->steps);
        free(edges[i]);
    }
    fclose(write_fptr);
}

float calc_pageRank(int v1, int v2, int N_VERTICES, float* prob_tran, float alpha, int itr){
    float transition_mat[N_VERTICES+1][N_VERTICES+1];
    for(int i=0;i<N_VERTICES+1;i++){
        for(int j=0;j<N_VERTICES+1;j++){
            transition_mat[i][j] = (float)prob_tran[idx(i, j, 1, N_VERTICES)];
        }
    }
    for(int i=1;i<N_VERTICES+1;i++){
        transition_mat[v2][i] *= (1-alpha);
        transition_mat[i][v1] = alpha;
    }
    float eig_vec[N_VERTICES+1];
    for(int i=0;i<N_VERTICES+1;i++){
        eig_vec[i] = 1.0;
    }
    while(itr>0){
        float xyz[N_VERTICES+1];
        float max = 0;
        for(int i=1;i<N_VERTICES+1;i++){
            xyz[i] = 0;
            for(int j=1;j<N_VERTICES+1;j++){
                xyz[i] += eig_vec[j]*transition_mat[i][j];
            }
            if(max < xyz[i]){
                max = xyz[i];
            }
        }
        for(int i=1;i<N_VERTICES+1;i++){
            eig_vec[i] = (float)xyz[i]/(float)max;            
        }
        itr--;
    }
    float sum = 0;
    for(int i=1;i<N_VERTICES+1;i++){
        sum += eig_vec[i];
    }
    for(int i=1;i<N_VERTICES+1;i++){
        eig_vec[i] /= sum;
    }
    return eig_vec[v2];    
}


void pageRank(int N_VERTICES, int N_EDGES, float* prob_tran, struct Node* adj[], float alpha, int itr, int K){
    int total_edges = (N_VERTICES*(N_VERTICES-1)/2) - N_EDGES;
    struct Edge* edges[total_edges];
    int cnt = 0;
    for(int i=1;i<N_VERTICES + 1;i++){
        for(int j=i+1;j<N_VERTICES + 1;j++){
            int flag = 0;
            struct Node* curr = adj[i];
            while(curr){
                if(curr->data == j){
                    flag = 1;
                    break;
                }
                curr = curr->next;
            }
            if(flag == 0){
                float score = calc_pageRank(i, j, N_VERTICES, prob_tran, alpha, itr);
                score += calc_pageRank(j, i, N_VERTICES, prob_tran, alpha, itr);
                score /= 2;
                edges[cnt] = newEdge(i, j, score);
                cnt++;
            }
        }
    }
    qsort(edges, total_edges, sizeof(struct Edge*), comparator_function);
    FILE *write_fptr = fopen("PageRank.txt", "w");
    for(int i=0;i<K;i++){
        fprintf(write_fptr, "%d %d %f\n", edges[i]->first, edges[i]->second, edges[i]->score);
        free(edges[i]);
    }
    fclose(write_fptr);

}

// Function to compute the adjacency list and adjacency matrix (required and used only for Katz score computation) of the graph
void compute_adjacency_list_and_matrix(int N_EDGES, int N_VERTICES, FILE* fptr, struct Node* adj[], long A[][N_VERTICES+1]){
    if(fptr == NULL){
        printf("Error opening file");
        exit(1);
    }
    for(int i=0;i<N_VERTICES+1;i++){
        adj[i] = malloc(sizeof(struct Node));
        adj[i]->data = i;
        adj[i]->next = NULL;
    }
    for(int i=0;i<N_VERTICES+1;i++){
        for(int j=0;j<N_VERTICES+1;j++){
            A[i][j] = 0;
        }
    }
    int v1, v2, w1;
    int c = 0;
    while(c < N_EDGES){
        fscanf(fptr, "%d %d %d", &v1, &v2, &w1);
        addEdgeUDG(v1, v2, adj);
        A[v1][v2] = 1;
        A[v2][v1] = 1;
        c++;
    }
    printGraph(N_VERTICES, adj);
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

// Function to count number of vertices in a graph
int count_vertices(int N_EDGES, FILE* fptr){
    if(fptr == NULL){
        printf("Error opening file");
        exit(1);
    }
    int v1, v2, w1;
    int c = 0;
    struct Node* vertex_list = newNode(-1);
    while(c < N_EDGES){
        fscanf(fptr, "%d %d %d", &v1, &v2, &w1);
        add_vertex(vertex_list, v1);
        add_vertex(vertex_list, v2);
        c++;
    }
    return length_ll(vertex_list->next);
}

int main(){
    
    // Input file name
    char file_name[] = "../contact-high-school-proj-graph.txt";
    
    // Find out how many edges and vertices are there in the graph
    FILE* fptr = fopen(file_name, "r");
    int N_EDGES = count_edges(fptr);
    fclose(fptr);
    fptr = fopen(file_name, "r");
    int N_VERTICES = count_vertices(N_EDGES, fptr);
    fclose(fptr);
    
    // Initialize and create adjacency list and matrix
    struct Node* adj[N_VERTICES+1];
    long A[N_VERTICES+1][N_VERTICES+1];

    // Compute adjacency list
    // Also compute adjacent matrix (required in Katz score calculation)(not used for any other part)
    fptr = fopen(file_name, "r");
    compute_adjacency_list_and_matrix(N_EDGES, N_VERTICES, fptr, adj, A);
    fclose(fptr);
    printf("Adjacency list printed. Check adjList.txt.\n");

    int K;
    printf("Enter Value for K -> ");
    scanf("%d", &K);

    // Compute Jaccard Scores
    jaccard(N_EDGES, N_VERTICES, adj, K);
    printf("Jaccard scores for top %d edges printed. Check Jaccard.txt.\n", K);

    // Compute Katz Scores
    long* path_count = malloc(sizeof(path_count) * (N_VERTICES+1) * (N_VERTICES+1) * 7);
    float beta = 0.1;
    katz(N_VERTICES, N_EDGES, adj, A, path_count, beta, K);
    free(path_count);
    printf("Katz scores for top %d edges printed. Check Katz.txt.\n", K);

    // Initialization for the next two parts
    int steps = 1000;
    float* prob_tran = malloc(sizeof(prob_tran) * (N_VERTICES+1) * (N_VERTICES+1) * steps);
    
    // Compute Commute time
    steps = 6;
    commuteTime(steps, N_VERTICES, N_EDGES, prob_tran, adj, K);
    printf("Commute time for top %d edges printed. Check HittingTime.txt.\n", K);

    // Compute Commute time Accurately
    steps = 50;
    commuteTimeAccurate(steps, N_VERTICES, N_EDGES, prob_tran, adj, K);
    // free(prob_tran);
    printf("Commute time accurately for top %d edges printed. Check HittingTimeAccurate.txt.\n", K);

    // Compute Page Rank
    float alpha = 0.2;
    int itr = 2; // Number of iterations of the power method
    pageRank(N_VERTICES, N_EDGES, prob_tran, adj, alpha, itr, K);
    free(prob_tran);
    printf("Rooted Page Rank for top %d edges printed. Check PageRank.txt.\n", K);

    for(int i=0;i<N_VERTICES+1;i++){
        free(adj[i]);
    }
    return 0;
}
