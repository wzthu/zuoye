# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 19:22:29 2018
@author: WeiZheng
"""

import numpy as np

class ValidGraph:
    def __init__(self,edges,nodes):
        self.nodeIdx = self.genNodeIdx(edges, nodes)
        self.nodeAttr = {}        
        self.adjMatrix = self.getAdjMatrix(self.nodeIdx,edges)
        
    def genNodeIdx(self,edges, nodes):
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
        for node in nodes:
            if not isinstance(node,str):
                raise Exception("node",node,"should be an string")
            if node not in nodeIdx:
                nodeIdx[node] = idx
                idx +=1
                    
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
    def __init__(self,graphsEdgeList, graphsNodeList):
        if len(graphsEdgeList) != len(graphsNodeList):
            raise Exception('the number of list in graphsEdgeList and graphsNodeList should be the same as the number of graphs')
        graphNumb = len(graphsEdgeList)
        self.graphs = []
        mergeEdges = []
        mergeNodes = []
        for i in range(graphNumb):   
            edges = graphsEdgeList[i]
            nodes = graphsNodeList[i]
            self.graphs.append(ValidGraph(edges, nodes))  
            mergeEdges.extend(edges)
            mergeNodes.extend(nodes)
      
        self.mergeGraph = ValidGraph(mergeEdges,mergeNodes) 
    
    def isConnect(self,upstream,downsteam,paramNum):
        if paramNum >=0 and paramNum < len(self.graphs):
            return self.graphs[paramNum].isConnect(upstream,downsteam)
        else:
            raise Exception("paramNum",paramNum,'is not available, max paramNum is',len(self.graphs)-1)


class GraphAll(GraphMng):
    def __init__(self,*args):
                #Smart-seq
        node1 = ['SRAToFastq',
                 'FastQC',
                 'FastqDump',
                 'Hisat2',
                 'Tophat',
                 'Star',
                 'Cufflinks',
                 'HTSeq',
                # 10x涉及的类的名称
                 'Quantification10x',
                 'PCA',
                # scATAC-seq
                 'FastQC',
                 'AdapterRemoval',
                 'Bowtie2',
                 'SamToBam',
                 'BamToSam',
                 'DuplicateRemoval'                 
                ]
        
        edge1 = [
                #Smart-seq
                ['FastqDump','Hisat2'],
                ['SRAToFastq','FastQC'],
                ['SRAToFastq','Tophat'],
                ['SRAToFastq','Star'],
                ['Tophat','Cufflinks'],
                ['Star','HTSeq'],
                #10x Genomeics
                ['Cellranger','Seurat'],
                #drop-seq
                ['FastqToBam','BamMerge'],
                ['BamMerge','TagBarcode'],
                ['TagBarcode','TagBarcode'],
                ['TagBarcode','TrimAdapter'],
                ['TrimAdapter','TrimPolyA'],
                ['TrimPolyA','BamToFastq'],
                ['BamToFastq','StarAlign'],
                ['StarAlign','StarBam'],
                ['StarBam','MergeBamAlign'],
                ['TrimPolyA','MergeBamAlign'],
                ['MergeBamAlign','TagGene'],
                ['TagGene','DetectError'],
                ['DetectError','DigitalExpression'],
                #ATAC-seq
                ['SRAToFastq','FastQC'],
                ['SRAToFastq','AdapterRemoval'],
                ['AdapterRemoval','Bowtie2'],
                ['Bowtie2','SamToBam'],
                ['SamToBam', 'BamSort'],
                ['BamSort', 'RmDuplicates'],
                ['RmDuplicates', 'BamToBed'],
                ['BamToBed', 'MergeToFrag'],
                ['BamToBed', 'RmChrOrMergeAllSample'],
                ['RmChrOrMergeAllSample', 'MergeToFrag'],
                ['RmChrOrMergeAllSample', 'BedSort']
                ]
        super(GraphAll, self).__init__([edge1],[node1])        
        
class GraphATACgl(GraphMng):
    def __init__(self,*args):
        edge1 = [
                ['UnzipAndMerge','AdapterRemoval'],
                ['UnzipAndMerge','FastQC'],
                ['AdapterRemoval','Bowtie2'],
                ['Bowtie2','SamToBam'],
                ['BamToSam','SamToBam'],
                ['SamToBam','BamToSam'],
                ]
        super(GraphATACgl, self).__init__([edge1],[[]])
        
                
        
        
        