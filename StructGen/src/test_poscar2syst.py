from StructGen import StructGen
import kwant
import numpy as np 
from glob import glob 

POSCARS = glob('./POSCARS/*') 

gen = StructGen()

for POS in POSCARS: 
    print(POS)
    gen.poscar2syst(POS)
    syst=gen.get_syst()
    #kwant.plot(syst)
