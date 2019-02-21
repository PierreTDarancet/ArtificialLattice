from math import floor
Zigzag = kwant.lattice.general([[1,0],[0,np.sqrt(3)/3]], #Lattice vectors
                                     [[1/6,np.sqrt(3)/2],[2/6,0],[4/6,0],[5/6,np.sqrt(3)/2]]) # Coordinates
def get_width(N=7): 
    if N < 2: 
        raise("N cannot be less than 2")
    else:
        return int(N/2)*Zigzag.prim_vecs[1][1]
    
def get_length(L=8): 
    if L < 2: 
        raise("L cannot be less than 1")
    else:
        return (L/4)*Zigzag.prim_vecs[0][0]

def delete_section(pos): 
    x,y = pos 
    a = Zigzag.prim_vecs[0][0]
    b = Zigzag.prim_vecs[1][1]
    n = floor(x/a)
    if 0< y < n*1.5*b: 
        return True 
    else: 
        return False
    
def make_1D_zigzag(N=7):
    syst = kwant.Builder(kwant.TranslationalSymmetry(Zigzag.prim_vecs[0]))
    syst[Zigzag.shape((lambda pos: pos[1] >0 and pos[1] <= get_width(N)),(0,0))] = 0
    syst[Zigzag.neighbors()] = -1
    return syst
