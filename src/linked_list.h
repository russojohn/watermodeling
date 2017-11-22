#ifndef LL_H
#define LL_H

typedef struct _node {
	int data;
	struct _node * next;
} node;

void print_list(node* head);

void add_list(node** head, int data);

void remove_list(node** head, int data);
	


#endif // !LL_H
