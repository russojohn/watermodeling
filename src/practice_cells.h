#ifndef CELLS_H
#define CELLS_H


typedef struct _listcell {
	int *HoC;                            // Head of Chain for linked list
	int *LinkedList;                     // linked list
	int NumberCells_x;                     // number of cells in one direction
	int NumberCells_y;
	int NumberCells_z;
	double CellSize_x;                     // cell edge size
	double CellSize_y;
	double CellSize_z;
} listcell;


listcell* getList(vector box_sides,double cutoff,int num_particles);
void freeList(listcell *l);
void resetList(listcell *l);
void updateList(listcell *l,const vector *pos,int num);
int getParticleEnergy(listcell *l,vector *pos,int num_particles,int label,vector box, double *energy);
void changeCell(listcell *l,const vector *oldpos,const vector *newpos,int num);

#endif
