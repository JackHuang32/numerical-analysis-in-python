import numpy as np
from numpy import*
from numpy import linalg
def current_pivot(M,sv,t,r):
    m,n = np.shape(M)
    max = 0
    idx = -1
    for i in range(m):
        if i not in r:
            if max < abs(M[i][t])/sv[i]:
                max = abs(M[i][t])/sv[i]
                idx = i
    r.append(idx)
    return idx
def find_last_row(m,r):
    for i in range(m):
        if i not in r:
            return i
def SPP(A,b):
    print('original system:\n',np.column_stack((A,b)))
    m_A,n_A = np.shape(A)
    sv = [max(abs(A[i])) for i in range(m_A)]
    augmented_matrix = np.column_stack((A,b))
    m,n = np.shape(augmented_matrix)
    L = np.zeros((m,m))
    solution = [-1 for i in range(m)]
    r = []   #for checking if the row is selected and record the permutation
    for t in range(m-1):#eliminate m-1 times 
        idx = current_pivot(augmented_matrix,sv,t,r)#find current pivot row
        
        for i in range(m):#eliminate all rows except for the current and above rows
            if i not in r:
                rate = augmented_matrix[i][t]/augmented_matrix[idx][t]
                for j in range(t,n):
                    augmented_matrix[i][j]-=(rate*augmented_matrix[idx][j])
                L[i][t]=rate
                
                
    last_row = find_last_row(m,r)#find the row haven't been selected
    r.append(last_row)
    Lower = copy(L)
    for i in range(m):
        Lower[i] = L[r[i]]
    Lower = np.add(Lower,identity(m))
    print('r:',r)
    print('current system: \n',augmented_matrix)
    print('upper triangular:\n')
    upper_triangular = copy(augmented_matrix)
    for i in range(m):
        upper_triangular[i] = augmented_matrix[r[i]]
    print(upper_triangular)
    print('lower triangular:\n',Lower)
    #start of the backward subsitution
    
    backward_col_idx = n_A
    for i in range(len(r)-1,-1,-1):#follow the record of elimination order backward
        current_row_sum = 0
        for j in range(backward_col_idx,n_A):
            current_row_sum+=augmented_matrix[r[i]][j]*solution[j]
        solution[backward_col_idx-1] = (augmented_matrix[r[i]][n-1]-current_row_sum)/augmented_matrix[r[i]][i]
        backward_col_idx -= 1
    print('solution:\n',solution)
A = np.array([[-9,11,-21,63,-252],[70,-69,141,-421,1684],[-575,575,-1149,3451,-13801],[3891,-3891,7782,-23345,93365],[1024,-1024,2048,-6144,24572]],dtype=double)
b = np.array([-356,2385,-19551,132274,34812],dtype=double)
SPP(A,b)