mat=[[1,2,3,0],
     [4,5,6,45],
     [7,8,9,83]]
def create(k,l,m,n,c):
    l1,l2,l3,l4=[],[],[],[]
    for i in range(k,n-1):
            l1.append(mat[k][i])   #creates a list from 1 to 3
            mat[k][i]=None
    for i in range(l,m-1):
            l2.append(mat[i][n-1]) #creates a list from 0 to 45
            mat[i][n-1]=None
    for i in range(n-1,k,-1):
            l3.append(mat[m-1][i]) #creates a list from 83 to 8
            mat[m-1][i]=None
    for i in range(m-1,l,-1):
            l4.append(mat[i][l])   #creates a list from 7 to 4
            mat[i][l]=None
    if(c==1):
        return(l1+l2+l3+l4)  #1 when clockwise spiral starts from top-left corner
    elif(c==2):
        return(l2+l3+l4+l1)  #2 when spiral starts from top-right corner
    elif(c==3):
        return(l3+l4+l2+l1)  #3 when spiral starts from bottom right corner
    elif(c==4):
        return(l4+l1+l2+l3)  #4 when spiral starts from bottom left corner
o=[]
for i in range(2): #2 is no of columns/2
    o+=([j for j in create(i,i,3,4,4)])
print([i for i in o if i is not None])    
    
#this way all cases can be considered, the key will tell from where to start
#for anti-clockwise i will use reverse
#also clockwise inwards will be reverse of anticlockwise outwards spiral, so inwards and outwards both can be considered