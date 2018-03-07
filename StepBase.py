# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 14:05:43 2018

@author: WeiZheng
"""

import os
from hashlib import md5
from GraphMng import GraphAll,GraphATACgl
import time
import subprocess




class Configure:   
    __config = {
        'threads' : 2,
        'genome' : None,
        'tmpdir' : os.path.abspath('./'),
        'refdir' : None,
        'pipeline': 'All',
        'regpipe':{
                'All':GraphAll(),
                'ATACgl':GraphATACgl()        
                },        
        }   
   
    def __init__(self,):
        raise Exception("can not be instanized")
    @classmethod
    def getConfig(cls,key):
        return cls.__config[key]
    
    @classmethod
    def setConfig(cls,key,val):
        if key == 'threads':
            Configure.setThreads(val)
        elif key == 'genome':
            Configure.setGenome(val)
        elif key == 'tmpdir':
            Configure.setTmpDir(val)    
        elif key == 'refdir':
            Configure.setRefDir(val)
        elif key == 'pipeline':
            Configure.setPipe(val) 
        elif key == 'regpipe':
            raise Exception('not configurable')
        else:
            cls.__config[key] = val
    
    @classmethod
    def setThreads(cls,val):
        cls.__config['threads']
    
    @classmethod
    def getThreads(cls,):
        # get the globle configure of threads size
        return cls.__config['threads']
    
    @classmethod
    def setGenome(cls,val):
        if cls.__config['refdir'] is None:
            raise Exception('refdir should be configure first, call Configure.setRefDir for configuration')
        ## bowtie2 index
        cls.__config['bt2Idx'] = os.path.join(cls.getRefDir(),val)
        suffix = ['.1.bt2','.2.bt2','.3.bt2','.4.bt2','.rev.1.bt2','.rev.2.bt2']
        cls.__config['bt2IdxFiles'] = [ cls.__config['bt2Idx'] + s for s in suffix ]
        
        
    @classmethod
    def getGenome(cls):
        return cls.__config['genome']
    
    @staticmethod
    def checkFolderPath(folderPath):
        if not os.path.isdir(os.path.abspath(folderPath)):
            raise Exception(folderPath,"is not an folder")
        if not os.path.exists(folderPath):
            raise Exception(folderPath,"is not exist")
        if not (os.access(folderPath,os.X_OK) and os.access(folderPath,os.W_OK)):
            raise Exception(folderPath,"is not accessable")
        return True
    
    @classmethod
    def setRefDir(cls,folderPath):
        Configure.checkFolderPath(folderPath)
        cls.__config['refdir'] = folderPath
    
    @classmethod
    def getRefDir(cls,):
        return cls.__config['refdir']    
    
    @classmethod
    def setTmpDir(cls,folderPath):
        Configure.checkFolderPath(folderPath)
        cls.__config['tmpdir'] = folderPath
        
    @classmethod
    def getTmpDir(cls,):
        return cls.__config['tmpdir'] 
    
    @classmethod
    def getTmpPath(cls,foldOrFileName):        
        if isinstance(foldOrFileName,list):
            result = []
            for name in foldOrFileName:
                result.append(os.path.join(cls.getTmpDir(),name))
            return result
        else:
            return os.path.join(cls.getTmpDir(),foldOrFileName)
           
    @classmethod
    def setPipe(cls,pipeName):
        if not pipeName in cls.__config['regpipe']:
            raise Exception(pipeName,':pipeName not found')
        cls.__config['pipeline'] = pipeName
        
    @classmethod
    def getPipe(cls,):
        return cls.__config['pipeline']
    
    @classmethod
    def addPipe(cls,pipeName,graphMngObj):
        cls.__config['regpipe'][pipeName] = graphMngObj
    
    @classmethod
    def getGraph(cls,pipeName = None):
        if pipeName is None:            
            return cls.__config['regpipe'][cls.getPipe()]
        else:
            return cls.__config['regpipe'][pipeName]
    
class Schedule:
    __schedule = []
    
    @classmethod
    def add(cls, stepObj):
        if isinstance(stepObj,Step):
            cls.__schedule.append(stepObj)
        else:
            raise Exception('only support schedule sub-classes')
    @classmethod
    def run(cls, ):
        for step in cls.__schedule:
            step.run()
    
    @classmethod
    def remove(cls,step):       
        if isinstance(step,Step):
            stepid = -1
            stepid = step.getStepID()
            for i in range(len(cls.__schedule)):
                if stepid == cls.__schedule[i].getStepID:
                    del(cls.__schedule[i])
        elif isinstance(step,int):
            del(cls.__schedule[step])
            
    @classmethod
    def getSchedule(cls,):
        return cls.__schedule
        
            
    
        

class StepBase:  
    stepObjCounter = 0
    @classmethod
    def regStepID(cls):
        objid = cls.stepObjCounter
        cls.stepObjCounter += 1
        return objid
    
    def __init__(self,cmdParam,**kwargs):
        self.inputs = {}    
        self.outputs = {}
        self.paramsIO = {}
        self.params = {}
        self.unsetParams = ''
        self.funGroup = None
        self.__isFinished = False
        self.__multiRun = False
        self.checkParamValid()
        self.__stepID = StepBase.regStepID()
        self.__inputSize = -1
        self.__upstreamSize = 1
        return 0
    
    def getStepID(self,):
        return self.__stepID
    
    def getStepFolderName(self,):
        return 'step_' + str(self.getStepID()).zfill(2) + '_' + self.__class__.__name__
    
    def getStepFolerPath(self,):
        return Configure.getTmpPath(self.getStepFolderName())
    
    def initIO(self,): 
        tmpdir = Configure.getTmpDir()
        if not os.path.exists(self.getStepFolerPath()):
            os.mkdir(self.getStepFolerPath())
        Configure.setTmpDir(self.getStepFolerPath())
        self.impInitIO()
        Configure.setTmpDir(tmpdir)
        self.checkRunnable()
        
    def impInitIO(self,):
        raise Exception("method: initIO must be overwrited")
        
    def __call__(self, *args): 
        notNoneNumb = 0
        for i in range(len(args)):
            if not args[i] is None:
                notNoneNumb += 1
            self.validUpstream(stepObj = args[i],paramNum = i)
        if notNoneNumb == 0:
            raise Exception('there should be one upstream Step object in parameter')
        if len(args) != self.__upstreamSize:
            raise Exception("only support", self.__upstreamSize,"upstream")
        self.call(*args)
        self.initIO()
        return self
      
    ## subclass implement
    ## object order can not be changed
    def call(self, *args):
        raise Exception("method: call must be overwrited")
        return 0
    
    def getFileNamePrefix(self,fileName):
        return os.path.split(fileName)[-1]
    
    def getMaxFileNamePrefix(self,file1,file2):
        file1 = self.getFileNamePrefix(file1)
        file2 = self.getFileNamePrefix(file2)
        len1 = len(file1)
        len2 = len(file2)        
        for i in range(min(len1,len2)):
            if len1[i] != len2[i]:
                break        
        if i == 0:
            return ''
        elif len1[i-1]=='.':
            return len1[0:i-1]
        else: 
            return len1[0:i]
    
    def validUpstream(self,stepObj,paramNum=0):
        if stepObj is None:
            return True
        if not isinstance(stepObj,Step):
            raise Exception("type error")
        if not Configure.getGraph().isConnect(stepObj.__class__.__name__,self.__class__.__name__,paramNum):
            raise Exception("")
        stepObj.checkFilePath(checkExist = False)
        return True
    
    def checkRunnable(self,addToSchedule = True):
        try:
            self.checkFilePath(checkExist=False)
            Schedule.add(self)
            print(self.getStepFolderName +'is ready and added to Schedule')
        except Exception:
            pass
        
        
    def getLogName(self,):
        return self.__class__.__name__ + '.' + self.getParaMD5code() + '.log'
    
    def getLogPath(self,):
        return os.path.join(self.getStepFolerPath(),self.getLogName())
    
    def _writeLogLines(self,strlines):
        if not os.path.exists(self.__logpath):
            raise Exception("can not write log when log file is not created")
        if not isinstance(strlines,list):
            strlines = [strlines]        
        logfile = open(self.__logpath,'a')
        logfile.writelines([ '||' + s +'\n' for s in strlines])
        logfile.close()
    
    def getCurTime(self,):
        return time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))
    
    
    def run(self,):
        self.checkInputFilePath()
        self.checkOutputFilePath(checkExist = False)
        if self.checkFinish():
            lines = ['=======================================',
            self.getCurTime()]
            self._writeLogLines(lines)
            print("Finished nothing to do")
            self.loadResult()
        else: 
            logpath = self.getLogPath()
            logFile = open(logpath,'w')
            logFile.close()
            self.__logpath = logpath
            lines = ['=======================================',
            self.getCurTime()]  
            self._writeLogLines(lines)
            tmpdir = Configure.getTmpDir()            
            Configure.setTmpDir(self.getStepFolerPath())
            if self.__multiRun :
                self._multiRun()
            else:
                if self.__inputSize == -1:
                    raise Exception('call self._setInputSize(your sample size) in impInitIO')
                for i in range(self.__inputSize):
                    self._singleRun(i)
            Configure.setTmpDir(tmpdir)
            if self.checkResult():
                self.setFinish()
        lines = [self.getCurTime(),
                '=======================================']
        self._writeLogLines(lines)
        return True
    
    ## subclass overwrite if not commandline
    def loadResult(self,):
        self.setFinish()
        return True
    
    
    ## one of the run method should be implement in subclass
    def _singleRun(self, i):
        raise Exception("method: _singleRun must be overwrited")
        
    def _multiRun(self,):
        raise Exception("method: _multiRun must be overwrited")
    

        
    

        
    def checkFinish(self,):  
        if self.__isFinished:
            return True
        logfilepath= Configure.getTmpPath(os.path.join(self.getStepFolderName(),self.getLogName()))
        if os.path.exists(logfilepath):
            self.__logpath = logfilepath
            try:
                self.checkOutputFilePath()
                return True
            except Exception:
                pass
            
        return False
    
    def setFinish(self,):        
        self.__isFinished = True
    
    def getFileNameAndSize(self,filePath,fileSize = True):
        namesizelist = []
        if not isinstance(filePath,list):
            filePath = [filePath]
        for filepath in filePath:
            filename = os.path.split(filepath)[-1]
            if fileSize: 
                filesize = os.path.getsize(filepath)
            else:
                filesize = ''
            namesizelist.append(filename + str(filesize))
        return namesizelist
    
    def getParaMD5code(self,):
        checklist = ['']
        checklist1 = []
        checklist2 = []
        checklist3 = []
        for key in self.inputs.keys():
            checklist1.extend(self.getFileNameAndSize(self.inputs[key]))        
        checklist1.sort()
        checklist.extend(checklist1) 
        for key in self.outputs.keys():
            checklist2.extend(self.getFileNameAndSize(self.outputs[key], fileSize = False))
        checklist2.sort()
        checklist.extend(checklist2)
        keys = list(self.params.keys())
        keys.sort()
        for key in keys:
            checklist3.append(self.params[key])
        checklist.extend(str(checklist3))
        unp = list(self.unsetParams.split())
        unp.sort()
        checklist.extend(unp)
        checklist = ''.join(checklist)
        return md5(checklist.encode('utf8')).hexdigest()[0:8]
        
        
            
        
    def checkFilePath(self, checkExist = True):
        self.checkOutputFilePath(checkExist)
        self.checkInputFilePath(checkExist)              
                
    
    def checkOutputFilePath(self, checkExist = True):
        for key in self.outputs.keys():
            outPaths = self.convertToList(self.outputs[key])
            self.checkFilePathList(outPaths,'output',key,checkExist)
    
    def checkInputFilePath(self, checkExist = True):
        for key in self.inputs.keys():
            inPaths = self.convertToList(self.inputs[key])
            self.checkFilePathList(inPaths,'input',key,checkExist)
    
    def convertToList(self,obj):
        if not isinstance(obj,list):
            return [obj]
        else:
            return obj
        
    
    def checkFilePathList(self,filePathList,iodict,key = None, checkExist = True):
        for filePath in filePathList:
            if filePath is None:
                raise Exception('File path of',iodict,key,' can not be None')
            if checkExist:
                if not os.path.exists(filePath):
                    raise Exception('File path of',iodict,key,' not found:',filePath)
                
    def getFileList(self, pathPrefix, filePath):
        fileList = self.getListInFile(filePath)
        for i in range(len(fileList)):
            fileList[i] =  os.path.join(pathPrefix,fileList[i])
        return fileList
    
    def getListInFile(self, filePath):
        listFile = open(filePath)
        fileList = listFile.readlines()
        if len(fileList) == 0:
            raise Exception('file',filePath,'can not be empty')
        return fileList
    
    def checkResult(self,):
        return True
    
    def checkParamValid(self,):
        return True
    
    def getInputs(self,):
        return list(self.inputs.keys())
           
    def getInput(self, inputName):
        return self.inputs[inputName]
    


    def setInput(self, inputName, inputValue):
        self.inputs[inputName] = inputValue

    
    
    def getOutputs(self,):
        return list(self.outputs.keys())
    
    def getOutput(self,outputName):
        return self.outputs[outputName]
    
    def getOutputList(self, outputName):
        if isinstance(self.outputs[outputName],list):
            return self.outputs[outputName]
        else:
            return [self.outputs[outputName]]
    
    def setOutput(self, outputName, outputValue):
        self.outputs[outputName] = outputValue
    
    def getUnsetParams(self,):
        return self.unsetParams
    
    def getParams(self,):
        return list(self.params.keys())
    
    def getParam(self, paramName):
        return self.params[paramName]
    
    def setParam(self, paramName, paramValue):
        """
        Set parameters (except for input or output parameters) in __init__()
        """
        self.params[paramName] = paramValue
        
    def getParamIOs(self,):
        """
        Get parameters (except for input or output parameters) setted in __init__()
        """
        return list(self.paramsIO.keys())
    
    def getParamIO(self, paramName):
        """
        Get input or output parameters set in __init__() or call()
        """
        return self.paramsIO[paramName]
    
    def setParamIO(self, paramName, paramValue):
        """
        Set input  or output parameters at __init__() and 
        Set input parameters from upstream at call()
        """
        if paramValue is None:            
            self.paramsIO[paramName] = None
        else:
            if isinstance(paramValue,list):
                self.paramsIO[paramName] = [os.path.abspath(s) for s in paramValue]
            else:
                self.paramsIO[paramName] = os.path.abspath(paramValue)
            
            
    
    def getConfigVal(key):
        return Configure.getConfig(key)
    
    def _setMultiRun(self,):
        self.__multiRun = True
        
    def _setInputSize(self,inputSize):
        self.__inputSize = inputSize
    
class Step(StepBase):
    def __init__(self,cmdParam,**kwargs):
        super(StepBase, self).__init__(cmdParam,**kwargs)
        
    def initInputWithFileList(self,inputName,inputPath,inputListFile):
        if inputListFile is None:           
            self.setInput(inputName, inputPath)
        elif inputPath is None:
            raise Exception(inputName,'should not be None')
        else:
            self.setInput(inputName,self.getFileList(inputPath,inputListFile))
    
    
    def setOutputDir1To1(self, outputName, outputDir, outputPrefix, outputSuffix, inputName):
        """
        Use the known input file paths to generate the output file paths.
        outputName(str): the key name of the file path list to be generate
        outputDir(str or None): the folder of the output file stored, None or string is OK.
                   if None, its value will set globally
        outputPrefix(str or None): the prefix name of output file name. 
                                   if it is None, the prefix of inputName will be use as the templete
        outputSuffix(str): the suffix name of output file name
        inputName(str): the input file paths to be referenced
        it will set paths list of outputDir/outputPrefix.*.outputSuffix like paths for outputName
        """
        inputList = self.getInput(inputName)
        if inputList is None:
            self.setOutput(outputName, None)
        else:
            if not isinstance(inputList,list):
                inputList = [inputList]
            outputList = [] 
            if outputPrefix is None:
                fileBasename = [os.path.basename(s) for s in inputList]
                prefixs = [s.split('.')[0] for s in fileBasename]
                for i in range(len(inputList)):
                    outputList.append(prefixs[i] + '.' + outputSuffix)
            else:
                for i in range(len(inputList)):
                    outputList.append(outputPrefix + '.' + str(i) + '.' + outputSuffix)
            if outputDir is None:
                self.setOutput(outputName,Configure.getTmpPath(outputList))
            else:        
                self.setOutput(outputName,[os.path.join(outputDir,s) for s in outputList])
                
    def setOutputDirNTo1(self, outputName, outputFilePath, fileNameForUnsetOutPath, inputName):
        """
        Use the known N input file paths to generate the one output file.
        outputName(str): the key name of the file path list to be generate
        outputFilePath(str or None): the path the output file, None or string is OK.
                   if None, its value will set globally
        fileNameForUnsetOutPath(str): the candidate name of output file name when outputFilePath is None. 
        inputName(str): the input file paths to be referenced        
        """
        inputList = self.getInput(inputName)
        if inputList is None:
            self.setOutput(outputName, None)
        else:
            if outputFilePath is None:
                if fileNameForUnsetOutPath is None:
                    raise Exception('fileNameForUnsetOutPath can not be None')
                else:
                    self.setOutput(outputName,Configure.getTmpPath(fileNameForUnsetOutPath))
            else:
                self.setOutput(outputFilePath)
            
    def callCmdline(self,cmdline,shell = False, stdoutToLog = True):
        if not shell:
            cmdline = ' '.join(cmdline)
        print(cmdline)
        try:
            if stdoutToLog:
                result = subprocess.run(cmdline,shell=True,check=True,stdout = subprocess.PIPE, stderr = subprocess.PIPE )
                self._writeLogLines(str(result.stdout))
                self._writeLogLines(str(result.stderr))
                return result
            else:
                result = subprocess.run(cmdline,shell=True,check=True) 
            
        except subprocess.CalledProcessError:
            raise 
            
    def getBoolParamCmd(self,cmdName,paramName):
        return cmdName if self.getParam(paramName) else ''
    
    def setInputDirOrFile(self, inputName, inputValue):
        """
        when input parameter is a list or single string,
        use this parameter to generate all of the file paths under the path or
        just use the single string path directly
        inputName: the key name of inputs
        inputValue: string of directory or file path, or list of file paths  
        """
        if inputValue is None:
            self.inputs[inputName] = None
        else:
            if isinstance(inputValue,list):
                self.inputs[inputName] = inputValue
            else:
                if os.path.isdir(inputValue):
                    filelist = os.listdir(inputValue)
                    filelist.sort()
                    suffix = [s.split('.')[-1] for s in filelist]
                    if len(set(suffix)) != 1:
                        raise Exception('the suffix of files under path:',inputName,'is not the same, check the file format under the directory')
                    self.inputs[inputName] = [os.path.join(inputValue,s) for s in filelist ]
                else:
                    self.inputs[inputName] = inputValue      
    
    def getInputList(self, inputName):
        """
        the wraper of getInput(). if the input is the single string,
        it will wrap it will list. Otherwise it will return the input as it is.
        """
        if isinstance(self.inputs[inputName],list):
            return self.inputs[inputName]
        else:
            return [self.inputs[inputName]]
        

        

        
