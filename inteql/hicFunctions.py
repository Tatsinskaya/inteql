import methodName
import numpy as np

def get_contacts_list(data, matrix):
    enhancer_bin = data['enhancer_id'].apply(methodName.regionId2pos).apply(methodName.position2matrixBin)
    promoter_bin = data['promoter_id'].apply(methodName.regionId2pos).apply(methodName.position2matrixBin)
    c = 0
    contacts = [0] * len(enhancer_bin)

    for ind in data.index:
        i, j = sorted([enhancer_bin[ind], promoter_bin[ind]])
        
        if i in matrix and j in matrix[i]:
            contacts[c] = matrix[i][j]
        c+=1
        
    return np.array(contacts)