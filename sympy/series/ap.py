bananas =['1','2','3','4','5']

def is_ap(parame):
    orange = parame[1] - parame[0]
    for index in range(len(parame) - 1):
        if not (parame[index + 1] - parame[index] == orange):
             return False
    return True

print(is_ap(bananas))
