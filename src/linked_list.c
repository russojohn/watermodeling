#include "linked_list.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void print_list(node* head) {
	node * curr = head;

	printf("printing list: ");
	while (curr != 0) {
		printf("%d\n ", curr->data);
		curr = curr->next;
	}
	printf("\n");
}

void add_list(node** head, int data){

	if (*head == 0){
		//init list
		*head = malloc(sizeof(node));

		(*head)->next = 0;
		(*head)->data = data;
	}
	else{
		//add to start
		node* temp = malloc(sizeof(node));

		temp->data = data;
		temp->next = *head;
		*head = temp;
	}
}

void remove_list(node** head, int data){
	//check if list is empty
	if (*head == 0){
		printf("List empty");
	}
	else if ((*head)->data == data){
		//remove from start
		node* temp_ptr;

		temp_ptr = (*head)->next;
		free(*head);
		*head = temp_ptr;
	}
	else{
		//remove from end/middle
		node *curr = (*head)->next, *prev = *head;

		while (curr != 0){
			if (curr->data == data){
				prev->next = curr->next;
				free(curr);
				return;
			}
			prev = curr;
			curr = curr->next;
		}

		printf("data does not belong to list\n");

	}

}

