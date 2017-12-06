import ovito
import math
from ovito.io import import_file
from ovito.modifiers import *
import subprocess
import sys

if len(sys.argv)!=2:
	print("%s [file]" % (sys.argv[0]))
	sys.exit()


pfile=open(sys.argv[1],"r")
line=pfile.readline()

box_x=float(line.split()[2])
box_y=float(line.split()[3])
box_z=float(line.split()[4])

pfile.close()

stringa="~/bin/convertxyz "+sys.argv[1]+" > xxx.xyz"
subprocess.Popen("~/bin/convertxyz "+sys.argv[1]+" > xxx.xyz",shell=True).wait()


node=import_file("xxx.xyz",columns =["Position.X", "Position.Y", "Position.Z"])

cell = node.source.cell

realcell = cell.matrix.copy()

realcell[2,0]=0
realcell[1,0]=0
realcell[0,0]=box_x

realcell[0,1]=0
realcell[2,1]=0
realcell[1,1]=box_y

realcell[0,2]=0
realcell[1,2]=0
realcell[2,2]=box_z

cell.matrix = realcell

node.add_to_scene()

node.modifiers.append(IdentifyDiamondModifier())
#node.modifiers.append(SelectExpressionModifier(expression="(StructureType!=1) && (StructureType!=4)"))
node.modifiers.append(SelectExpressionModifier(expression="StructureType==0"))
node.modifiers.append(DeleteSelectedParticlesModifier())

legami=CreateBondsModifier()
legami.cutoff=1.4
legami.bonds_display.width=0.12
node.modifiers.append(legami)


#pos = node.source.particle_properties.position      # ParticleProperty storing the positions
#pos.display.shape = ParticleDisplay.Shape.Square


