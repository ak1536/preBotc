'''
Created on Jul 14, 2016

@author: andrewkennedy
'''

def findClosest(vec, val):
    currentBest = 0
    lowness = lambda x : abs(x - val)
    bestValue = lambda : vec[currentBest]
    for i in range(len(vec)):
        if lowness(vec[i]) < lowness(bestValue()):
            currentBest = i
    return currentBest
    
if __name__ == '__main__':
    print findClosest(range(30), 12.7)

   

