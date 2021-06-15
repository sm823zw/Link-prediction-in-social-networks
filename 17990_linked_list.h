#include <stdio.h>
#include <stdlib.h>

// Data structure to store a node of a linked list
struct Node{
    int data;
    struct Node* next;
};

// Function to create a new node
struct Node* newNode(int data){
    struct Node* node = malloc(sizeof(struct Node));
    node->data = data;
    node->next = NULL;
    return node;
}

// Function to print a linked list
void print_list(struct Node* a){
    struct Node* y = a;
    while(y){
        printf("%d ", y->data);
        y = y->next;
    }
    printf("\n");
}

// Function to check if an element is present in a linked list or not
int member(struct Node* a, int val){
    struct Node* curr = a;
    while(curr){
        if(curr->data == val){
            return 1;
        }
        curr = curr->next;
    }
    return 0;
}

// Function to count neighbours
int count_neighbours(struct Node* a){
    struct Node* curr = a->next;
    int c = 0;
    while(curr){
        c++;
        curr = curr->next;
    }
    return c;
}

// Function to return the cardinality of union of two lists
int set_union(struct Node* a, struct Node* b){
    int c = 0;
    struct Node* curr = a;
    while(curr){
        c++;
        curr = curr->next;
    }
    curr = b;
    while(curr){
        if(!member(a, curr->data)){
            c++;
        }
        curr = curr->next;
    }
    return c;
}

// Function to return the cardinality of intersection of two lists
int set_intersection(struct Node* a, struct Node* b){
    b = b;
    int c = 0;
    struct Node* curr = a;
    while(curr){
        if(member(b, curr->data)){
            c++;
        }
        curr = curr->next;
    }
    return c;
}

// Function that returns length of a list
int length_ll(struct Node* vertex_list){
    int len = 0;
    struct Node* curr = vertex_list;
    while(curr){
        len++;
        curr = curr->next;
    }
    return len;
}

// Function that adds an element to the list only if not present
struct Node* add_vertex(struct Node* vertex_list, int v){
    struct Node* curr = vertex_list;
    while(curr->next){
        if(curr->data == v){
            return vertex_list;
        }
        curr = curr->next;
    }
    if(curr->data == v){
        return vertex_list;
    }
    curr->next = newNode(v);
    return vertex_list;
}
