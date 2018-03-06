# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 19:22:29 2018

@author: WeiZheng
"""

import numpy as np

class ValidGraph:
    def __init__(self,edges):
        self.nodeIdx = self.genNodeIdx(edges)
        self.nodeAttr = {}        
        self.adjMatrix = self.getAdjMatrix(self.nodeIdx,edges)
        
    def genNodeIdx(self,edges):
        nodeIdx = {}
        idx = 0
        for i in range(len(edges)):
            if len(edges[i]) != 2:
                raise Exception("edges",i,"do not contain 2 node")
            for j in range(2):
                if not isinstance(edges[i][j],str):
                    raise Exception("node",edges[i][j],"should be an string")
                if not edges[i][j] in nodeIdx:
                    nodeIdx[edges[i][j]] = idx
                    idx += 1
                    
        return nodeIdx
    
    def getAdjMatrix(self,nodeIdx,edges):
        nodeNum = len(nodeIdx)
        adjMatrix = np.zeros((nodeNum,nodeNum),dtype=bool)
        for i in range(len(edges)):
            adjMatrix[nodeIdx[edges[i][0]],nodeIdx[edges[i][1]]] = True
        return adjMatrix
        
    def setNodeAttrs(self,attrType,aDict):
        attrDict = self.nodeIdx.copy()
        for key in aDict.keys():
            attrDict[key] = aDict[key]
        self.nodeAttr[attrType] =  attrDict
    
    def isConnect(self,fromNode, toNode): 
        return self.adjMatrix[self.getNodeIdx(fromNode),self.getNodeIdx(toNode)]
    
    def getNodeIdx(self,nodeStr):
        return self.nodeIdx[nodeStr] 
    
    def getNodeNum(self,):
        return self.adjMatrix.shape[0]
     


class GraphMng:
    def __init__(self,graphsEdgeList):
        self.graphs = []
        mergeEdges = []
        for edges in graphsEdgeList:
            self.graphs.append(ValidGraph(edges))  
            mergeEdges.extend(edges)
        
        self.mergeGraph = ValidGraph(mergeEdges)  
        
    
    def isConnect(self,upstream,downsteam,paramNum):
        if paramNum >=0 and paramNum < len(self.graphs):
            return self.graphs[paramNum].isConnect(upstream,downsteam)
        else:
            raise Exception("paramNum",paramNum,'is not available, max paramNum is',len(self.graphs)-1)


class GraphAll(GraphMng):
    def __init__(self,*args):
        graph1 = [
                ['UnzipAndMerge','AdapterRemoval'],
                ['UnzipAndMerge','FastQC'],
                ['AdapterRemoval','Bowtie2'],
                ['Bowtie2','SamToBam'],
                ['BamToSam','SamToBam'],
                ['SamToBam','BamToSam'],
                ]
        super(GraphAll, self).__init__([graph1])        
        
class GraphATACgl(GraphMng):
    def __init__(self,*args):
        graph1 = [
                ['UnzipAndMerge','AdapterRemoval'],                
                ['UnzipAndMerge','FastQC'],
                ['AdapterRemoval','Bowtie2']
                ]
        super(GraphATACgl, self).__init__([graph1])
        
                
        
        
        
        
        
        
            
        