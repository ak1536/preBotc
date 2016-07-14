'''
Created on Jul 14, 2016

@author: andrewkennedy
'''

def findClosest(vec, val, start = 0):
    import matplotlib.pyplot as plt
    plt.plot(vec)
    plt.title('target = %s' % val)
    plt.show()
    for i in range(start, len(vec)):
        if abs(vec[i] - val) < abs(vec[i + 1] - val):
            return i
        else:
            vec[i] = vec[i + 1]
        
    pass
    
if __name__ == '__main__':
    print findClosest(range(30), 12.5)

