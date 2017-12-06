#!/usr/bin/env python

import sys
from math import *
import sys
import re
import copy
import numpy as np
import string

class Tetrahedron:
	# the default one is a top left tetrahedron
	
	def __init__(self,pos=[0.,0.,0.]):
		self.pos=np.array(pos)
		
		
		self.bonds=np.zeros((4,3))
		self.bonds[0]=np.array([0,0,1])  # up down direction
		self.bonds[1]=np.array([-sqrt(8./9.),0,-1./3.]) # left -right direction
		self.bonds[2]=np.array([sqrt(2./9.), sqrt(2./3.), -1./3.])
		self.bonds[3]=np.array([sqrt(2./9.),-sqrt(2./3.), -1./3.])
		
	def flip_horizontally(self):
		self.bonds[:,0]*=-1
	def flip_vertically(self):
		self.bonds[:,2]*=-1
	def shift(self,shift):
		self.pos+=shift
	
	def scale(self,factor):
		self.pos*=factor
	
	def print_pos(self):
		print ("%lf %lf %lf" % (self.pos[0],self.pos[1],self.pos[2]))


if len(sys.argv)!=5:
	print sys.argv[0]+" [stacking: ex. c3h2ch] [cells x] [cells y] [density]"
	sys.exit(1)

layer_string=sys.argv[1]
cells_x=string.atoi(sys.argv[2])
cells_y=string.atoi(sys.argv[3])
density=string.atof(sys.argv[4])


# BUILD THE STACKING SEQUENCE /////////////////////////////////////////////////////////
p = re.compile("[a-z]")

s_pos=[]
s_type=[]

for m in p.finditer(layer_string):
	s_pos.append(m.start())
	s_type.append(m.group())
	
stacking=[]

for i,(p,t) in enumerate(zip(s_pos,s_type)):
	
	if (i>0):
		prev_p+=1
		
		if (p!=prev_p):
			n_layers=int(layer_string[prev_p:p])
		
			for j in range(n_layers-1):
				stacking.append(prev_t)
			
	
	stacking.append(t)
	
	prev_p=p
	prev_t=t


p=len(layer_string)
prev_p+=1
if (prev_p!=p):
	n_layers=int(layer_string[prev_p:p])
	for j in range(n_layers-1):
		stacking.append(prev_t)

# ////////////////////////////////////////////////////


n_particles=0

# build the first layer

layer=[]
for i in range(cells_x):
	for j in range(cells_y):
		
		origin=[2*sqrt(2.)*i,2*sqrt(2./3.)*j,0]
		
		t1=Tetrahedron()
		t1.shift(origin)
		layer.append(copy.deepcopy(t1))
		
		t2=Tetrahedron()
		t2.shift(t1.pos+t1.bonds[3])
		t2.flip_horizontally()
		t2.flip_vertically()
		layer.append(copy.deepcopy(t2))
		
		t3=Tetrahedron()
		t3.shift(t2.pos+t2.bonds[1])
		layer.append(copy.deepcopy(t3))
		
		t4=Tetrahedron()
		t4.shift(t3.pos+t3.bonds[2])
		t4.flip_horizontally()
		t4.flip_vertically()
		layer.append(copy.deepcopy(t4))
		
		n_particles+=4
		

layers=[]

layers.append(layer)

for s in stacking:
	
	new_layer=[]
	
	if s=='h':
		for t in layers[-1]:
			if t.bonds[0][2]>0:
				
				t1=copy.deepcopy(t)
				t1.shift(t.bonds[0])
				t1.flip_vertically()
				new_layer.append(t1)
				
				t2=copy.deepcopy(t1)
				t2.shift(t1.bonds[1])
				t2.flip_horizontally()
				t2.flip_vertically()
				new_layer.append(t2)
				
				n_particles+=2
				
		
		
	if s=='c':
		for t in layers[-1]:
			if t.bonds[0][2]>0:
				t1=copy.deepcopy(t)
				t1.shift(t.bonds[0])
				t1.flip_vertically()
				t1.flip_horizontally()
				new_layer.append(t1)
				
				t2=copy.deepcopy(t1)
				t2.shift(t1.bonds[1])
				t2.flip_vertically()
				t2.flip_horizontally()
				new_layer.append(t2)
				
				n_particles+=2
		
	layers.append(new_layer)
			
			
volume=2*sqrt(2.)*cells_x*2*sqrt(2./3.)*cells_y*(4./3.)*(1.+len(stacking))

scale_factor=pow(float(n_particles)/volume/density,1./3.)


print "0 %d %lf %lf %lf" % (n_particles,scale_factor*2*sqrt(2.)*cells_x,scale_factor*2*sqrt(2./3.)*cells_y,scale_factor*4./3.*(1+len(stacking)))

for l in layers:
	for t in l:
		t.scale(scale_factor)
		t.print_pos()
