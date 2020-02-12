message = [[1,2,3,0],
          [4,5,6,45],
          [7,8,9,83]]
def create(k,l,m,n,c):
    l1,l2,l3,l4=[],[],[],[]
    for i in range(k,n-1):
            l1.append(message[k][i])   #creates a list from 1 to 3
            message[k][i]=None
    for i in range(l,m-1):
            l2.append(message[i][n-1]) #creates a list from 0 to 45
            message[i][n-1]=None
    for i in range(n-1,k,-1):
            l3.append(message[m-1][i]) #creates a list from 83 to 8
            message[m-1][i]=None
    for i in range(m-1,l,-1):
            l4.append(message[i][l])   #creates a list from 7 to 4
            message[i][l]=None
    if(c==1):
        return l1+l2+l3+l4  #1 when clockwise spiral starts from top-left corner
    elif(c==2):
        return l2+l3+l4+l1  #2 when spiral starts from top-right corner
    elif(c==3):
        return l3+l4+l2+l1  #3 when spiral starts from bottom right corner
    elif(c==4):
        return l4+l1+l2+l3  #4 when spiral starts from bottom left corner
'''1 when clockwise spiral starts from top-left corner
   2 when spiral starts from top-right corner
   3 when spiral starts from bottom right corner
   4 when spiral starts from bottom left corner'''


# will pass 1,2,3,4 according to route 
def encipher_route(message):
	i = int(0)
	j = int(0)
	temp = []
	while i<len(message) and j<len(message[0]):
		for k in create(i, j, len(message), len(message[0]), 1):
			if k is not None:
				temp.append(k)
		i = i+1
		j = j+1
	print(temp)

encipher_route(message)
