#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#include "vector.h"
#include "practice_energy.h"
#include "practice_cells.h"

static int module(int n,int mo)
{
	n=n%mo;
	
	while (n<0)
	{
		n+=mo;
	}
	
	return n;
}


listcell* getList(vector box_sides,double cutoff,int num_particles)
{
	listcell *l=malloc(sizeof(listcell));
	
	l->NumberCells_x=(int)(box_sides.x/cutoff);
	l->NumberCells_y=(int)(box_sides.y/cutoff);
	l->NumberCells_z=(int)(box_sides.z/cutoff);
	
	l->CellSize_x=box_sides.x/(double)l->NumberCells_x;
	l->CellSize_y=box_sides.y/(double)l->NumberCells_y;
	l->CellSize_z=box_sides.z/(double)l->NumberCells_z;
	
	l->HoC=(int*)calloc(l->NumberCells_x*l->NumberCells_y*l->NumberCells_z,sizeof(int));
	l->LinkedList=(int*)calloc(num_particles,sizeof(int));
		
	int i;
	int nnn=l->NumberCells_x*l->NumberCells_y*l->NumberCells_z;
	for (i=0;i<nnn;i++)
		(l->HoC)[i]=-1;
	
	return l;
}

void freeList(listcell *l)
{
	free(l->HoC);
	free(l->LinkedList);
	free(l);
}

void resetList(listcell *l)
{
	int i;
	int nnn=l->NumberCells_x*l->NumberCells_y*l->NumberCells_z;
	
	for (i=0;i<nnn;i++)
		(l->HoC)[i]=-1;
}

void updateList(listcell *l,const vector *pos,int num)
{
	int i;
	int ncell;                       // cell number
	int posx,posy,posz;              // cell coordinates
	int nnn;                         // total number of cells
	
	
	nnn=l->NumberCells_x*l->NumberCells_y*l->NumberCells_z;
	
	// HoC initialization
	for (i=0;i<nnn;i++)
		(l->HoC)[i]=-1;
	
	// colloids loop
	for (i=0;i<num;i++)
	{
		posx=(int)floor(pos[i].x/l->CellSize_x);
		posy=(int)floor(pos[i].y/l->CellSize_y);
		posz=(int)floor(pos[i].z/l->CellSize_z);
		
		posx=module(posx,l->NumberCells_x);
		posy=module(posy,l->NumberCells_y);
		posz=module(posz,l->NumberCells_z);
		
		ncell=posx+posy*l->NumberCells_x+posz*(l->NumberCells_x*l->NumberCells_y);
		
		(l->LinkedList)[i]=(l->HoC)[ncell];
		(l->HoC)[ncell]=i;
	}
	
}


int getParticleEnergy(listcell *l,vector *pos,int num_particles,int label,vector box, double *energy)
{
	int j;
	int nx,ny,nz,nn;               // Number of cells on surface and in box
	int neighbour_cell;                 // Neighbour cell
	int ix,iy,iz;                  // present cell coordinates
	int dx,dy,dz;
	int xv[3];                     // xlist values
	int yv[3];
	int zv[3];
	
	nx=l->NumberCells_x;
	ny=l->NumberCells_y;
	nz=l->NumberCells_z;
	nn=l->NumberCells_x*l->NumberCells_y;
	
	pos+=label;
	
	// reconstruct cartesian coordinates
	ix=(int)floor(pos->x/l->CellSize_x);
	iy=(int)floor(pos->y/l->CellSize_y);
	iz=(int)floor(pos->z/l->CellSize_z);
	
	pos-=label;
	
	// the next step shoudn't be necessary
	xv[0]=module(ix,nx);
	yv[0]=module(iy,ny);
	zv[0]=module(iz,nz);
	
	
	// cell number
	//int ncell=xv[0]+yv[0]*l->NumberCells_x+zv[0]*(l->NumberCells_x*l->NumberCells_y);
	
	
	xv[1]=module(ix+1,nx);
	yv[1]=module(iy+1,ny);
	zv[1]=module(iz+1,nz);
	
	
	xv[2]=module(ix-1,nx);
	yv[2]=module(iy-1,ny);
	zv[2]=module(iz-1,nz);
	
	
	// interactions with particles in the same cell
	// and in neighbouring cells
	
	
	double virial=0;
	int overlap=0;
	double e;
	int neighbours=0;
	
	for (dx=0;dx<3;dx++)
	{
		for (dy=0;dy<3;dy++)
		{
			for (dz=0;dz<3;dz++)
			{
				neighbour_cell=xv[dx]+yv[dy]*nx+zv[dz]*nn;
				
				j=(l->HoC)[neighbour_cell];
				
				while (j!=-1)
				{
					if (j!=label)
					{
						//computeparticleEnergy2Particles(label,j,pos,num_particles,box,energy,&virial);
						overlap=computeNumNeighbours(label,j,pos,num_particles,box,&e,&virial);
						
						if (overlap)
							return 1;
						else
							neighbours+=(int)e;
					}
					
					
					j=(l->LinkedList)[j];
				}
			}
		}
	}
	
	if (neighbours==4)
		*energy=-1;
	else
		*energy=0;
	
	return 0;
}



void changeCell(listcell *l,const vector *oldpos,const vector *newpos,int num)
{
	int ncell_old,ncell_new;         // cell number
	int posx,posy,posz;              // cell coordinates
	
	
	posx=(int)floor(oldpos->x/l->CellSize_x);
	posy=(int)floor(oldpos->y/l->CellSize_y);
	posz=(int)floor(oldpos->z/l->CellSize_z);
	
	posx=module(posx,l->NumberCells_x);
	posy=module(posy,l->NumberCells_y);
	posz=module(posz,l->NumberCells_z);
	
	ncell_old=posx+posy*l->NumberCells_x+posz*(l->NumberCells_x*l->NumberCells_y);
	
	posx=(int)floor(newpos->x/l->CellSize_x);
	posy=(int)floor(newpos->y/l->CellSize_y);
	posz=(int)floor(newpos->z/l->CellSize_z);
	
	posx=module(posx,l->NumberCells_x);
	posy=module(posy,l->NumberCells_y);
	posz=module(posz,l->NumberCells_z);
	
	ncell_new=posx+posy*l->NumberCells_x+posz*(l->NumberCells_x*l->NumberCells_y);
	
	if (ncell_old==ncell_new)
		return;
	
	int *old_item,*new_item;
	
	// delete the old position
	int j=(l->HoC)[ncell_old];
	
	old_item=l->HoC+ncell_old;
	
	while (j!=-1)
	{
		new_item=l->LinkedList+j;
		
		if (j==num)
		{
			*old_item=*new_item;
			break;
		}
		
		j=(l->LinkedList)[j];
		
		old_item=new_item;
	}
	
	// you should never be here because there must be a cancellation
	// in the list. This means that there are some problems in the list
	
	// add the new position
	(l->LinkedList)[num]=(l->HoC)[ncell_new];
	(l->HoC)[ncell_new]=num;
}

