import ovito
import math
from ovito.io import import_file
from ovito.modifiers import *
import subprocess
import sys

if len(sys.argv)!=2:
	print("%s [file]" % (sys.argv[0]))
	sys.exit()

stringa="~/bin/convertxyz "+sys.argv[1]+" > xxx.xyz"
subprocess.Popen("~/bin/convertxyz "+sys.argv[1]+" > xxx.xyz",shell=True).wait()


node=import_file("xxx.xyz",columns =["Position.X", "Position.Y", "Position.Z"])
node.add_to_scene()

node.modifiers.append(IdentifyDiamondModifier())
node.modifiers.append(SelectExpressionModifier(expression="(StructureType!=1) && (StructureType!=4)"))
node.modifiers.append(DeleteSelectedParticlesModifier())

legami=CreateBondsModifier()
legami.cutoff=1.4
legami.bonds_display.width=0.12
node.modifiers.append(legami)


#pos = node.source.particle_properties.position      # ParticleProperty storing the positions
#pos.display.shape = ParticleDisplay.Shape.Square


