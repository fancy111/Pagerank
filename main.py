import numpy as np
import pickle
import math
import os

fileIn = 'data/demo'
fileOut = 'data/demoSorted'

def getNodes():
    fin = open(fileIn, 'r')
    nodes = []
    for line in fin:
        words = line.split()
        source = int(words[0])
        target = int(words[1])
        if source not in nodes:
            nodes.append(source)
        if target not in nodes:
            nodes.append(target)
    fin.close()
    nodes.sort()
    print 'Nodes number:', len(nodes)
    return nodes

def mapNodes(nodes):
    nodesNum = len(nodes)
    src = nodes
    tar = [i for i in range(nodesNum)]
    res = dict(zip(src, tar))
    resRev = dict(zip(tar, src))
    print 'Map is created.'
    pickle.dump(res, open('map', 'w'))
    pickle.dump(resRev, open('mapRev', 'w'))

def sortData():
    mMap = pickle.load(open('map', 'r'))
    fin = open(fileIn, 'r')
    edges = []
    for line in fin:
        words = line.split()
        source = int(words[0])
        target = int(words[1])
        if [source, target] not in edges:
            edges.append([source, target])
    fin.close()
    edges = sorted(edges, key=lambda edge: edge[0])
    pickle.dump(edges, open('edges', 'w'))
    # edges = pickle.load(open('edges', 'r'))
    print 'Data is sorted.'

    fout = open(fileOut, 'w')
    for i in range(len(edges)):
        fout.write(str(mMap[edges[i][0]]) + '\t' + str(mMap[edges[i][1]]) + '\n')
    fout.close()

def getMatrix(dim, beta):
    A = np.zeros((dim, dim))
    fin = open(fileOut, 'r') # sorted
    srcOld = 0 # look at the start of the file
    src = srcOld
    outNodes = []
    for line in fin:
        words = line.split()
        src = int(words[0])
        tar = int(words[1])
        if src == srcOld:
            outNodes.append(tar)
        else:
            # print src, srcOld
            factor = 1.0 / len(outNodes)
            for i in outNodes:
                A[i][srcOld] = factor
            srcOld = src
            outNodes = []
            outNodes.append(tar)
    factor = 1.0 / len(outNodes)
    for i in outNodes:
        A[i][src] = factor
    fin.close()

    B = np.ones((dim, dim)) / dim
    A = beta * A + (1 - beta) * B
    pickle.dump(A, open('matrixA', 'w'))
    return A

def basic(nodesNum, A):
    r = np.array([1.0 / nodesNum for _ in range(nodesNum)])
    r = np.transpose(r)
    while True:
        rn = np.dot(A, r)
        # print rn
        e = np.linalg.norm((rn - r), ord=1)
        if e < 1e-6:
            return rn
        r = rn

def getSparseMatrix():
    fin = open(fileOut, 'r')  # sorted
    srcOld = 0  # look at the start of the file
    src = srcOld
    outNodes = []
    for line in fin:
        words = line.split()
        src = int(words[0])
        tar = int(words[1])
        if src == srcOld:
            outNodes.append(tar)
        else:
            tmp = []
            tmp.append(len(outNodes))
            for i in outNodes:
                tmp.append(i)
            pickle.dump(tmp, open('sparse/sparse_%d' % srcOld, 'w'))
            srcOld = src
            outNodes = []
            outNodes.append(tar)
    tmp = []
    tmp.append(len(outNodes))
    for i in outNodes:
        tmp.append(i)
    pickle.dump(tmp, open('sparse/sparse_%d' % srcOld, 'w'))
    fin.close()

def updateStep(nodesNum, beta, isFirstTime):
    if isFirstTime:
        r = np.array([(1.0 - beta) / nodesNum for _ in range(nodesNum)])
        f = open('rold/rold', 'w')
        pickle.dump(r, f)
        f.close()

    while True:
        r = pickle.load(open('rold/rold', 'r'))
        rn = np.array([(1.0 - beta) / nodesNum for _ in range(nodesNum)])
        for i in range(nodesNum):
            line = pickle.load(open('sparse/sparse_%d' % i, 'r'))
            di = line[0]
            destList = [nodes for nodes in line[1:]]
            for j in destList:
                rn[j] += beta * r[i] / di
        f = open('rold/rold', 'w')
        pickle.dump(rn, f)
        f.close()
        e = np.linalg.norm((rn - r), ord=1)
        if e < 1e-6:
            return rn

def getBlockStripeMatrix(nodesNum, basketSize):
    basket = int(math.ceil(float(nodesNum) / basketSize))
    fin = open(fileOut, 'r')  # sorted
    srcOld = 0  # look at the start of the file
    src = srcOld
    outNodes = []
    for line in fin:
        words = line.split()
        src = int(words[0])
        tar = int(words[1])
        if src == srcOld:
            outNodes.append(tar)
        else:
            tmp = [[len(outNodes)] for _ in range(basket)]
            for i in outNodes:
                tmp[i / basketSize].append(i)
            for i in range(basket):
                if len(tmp[i]) > 1:
                    pickle.dump(tmp[i], open('strip/strip_%d_%d' % (srcOld, i), 'w'))
                    print tmp[i]
            srcOld = src
            outNodes = []
            outNodes.append(tar)
    tmp = [[len(outNodes)] for i in range(basket)]
    for i in outNodes:
        tmp[i / basketSize].append(i)
    for i in range(basket):
        if len(tmp[i]) > 1:
            pickle.dump(tmp[i], open('strip/strip_%d_%d' % (srcOld, i), 'w'))
            print tmp[i]
    fin.close()

def blockStripeUpdate(nodesNum, basketSize, beta, isFirstTime):
    basket = int(math.ceil(float(nodesNum) / basketSize))

    if isFirstTime:
        for i in range(basket):
            r = np.array([(1.0 - beta) / nodesNum for _ in range(basketSize)])
            f = open('rold/rold_%d' % i, 'w')
            pickle.dump(r, f)
            f.close()

    while True:
        e = 0
        for i in range(basket):
            r = pickle.load(open('rold/rold_%d' % i, 'r'))
            rn = np.array([(1.0 - beta) / nodesNum for _ in range(basketSize)])
            for j in range(nodesNum):
                src = j
                if not os.path.exists('strip/strip_%d_%d' % (src, i)):
                    continue
                line = pickle.load(open('strip/strip_%d_%d' % (src, i), 'r'))
                di = line[0]
                destList = [nodes for nodes in line[1:]]
                for k in destList:
                    dest = k % basketSize
                    rn[dest] += beta * r[src % basketSize] / di
            f = open('rold/rnew_%d' % i, 'w')
            pickle.dump(rn, f)
            f.close()

            e += np.linalg.norm((rn - r), ord=1)
        if e < 1e-6:
            # print result
            x = []
            for i in range(basket):
                r = pickle.load(open('rold/rnew_%d' % i, 'r'))
                for i in r:
                    x.append(i)
                    x = x[:nodesNum]
            print x
            return

        for i in range(basket):
            rn = pickle.load(open('rold/rnew_%d' % i, 'r'))
            pickle.dump(rn, open('rold/rold_%d' % i, 'w'))

if __name__ == '__main__':
    nodes = getNodes()
    nodesNum = len(nodes)
    #
    # mapNodes(nodes)
    #
    # sortData()
    #
    A = getMatrix(nodesNum, 0.85)
    print basic(nodesNum, A)

    # getSparseMatrix()

    print updateStep(nodesNum, 0.85, True)

    basketSize = 2
    getBlockStripeMatrix(nodesNum, basketSize)
    blockStripeUpdate(nodesNum, basketSize, 0.85, True)

