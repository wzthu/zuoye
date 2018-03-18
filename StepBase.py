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
        'docker':True,
        'dockerPath':'/data',
        'dockerVersion':{
                'V1':'hca:latest',
                'V2':'hca:py2',
                #'V4':'hca:v4',
               # 'V':'hca:latest'
                },
        'identity':None,
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
        cls.__config['threads'] = val
    
    @classmethod
    def getThreads(cls,):
        # get the globle configure of threads size
        return cls.__config['threads']
    
    @classmethod
    def setDockerPath(cls,val):
        cls.__config['dockerPath'] = val
    
    @classmethod
    def getDockerPath(cls,):
        # get the globle configure of threads size
        return cls.__config['dockerPath']
    
    @classmethod
    def enableDocker(cls,val = True):
        if val:
            cls.getIdentity()
        cls.__config['docker'] = val
        
    @classmethod
    def isDocker(cls,):
        return cls.__config['docker']  
    
    @classmethod
    def setDockerVersion(cls,version,name):
        cls.__config['dockerVersion'][version] = name
    
    @classmethod 
    def getDockerVersions(cls):
        return cls.__config['dockerVersion'].keys()
    
    @classmethod 
    def getDockerVersion(cls,version):
        return cls.__config['dockerVersion'][version]
    
    @classmethod 
    def setIdentity(cls,identity):
        cls.__config['identity'] = identity
    
    @classmethod 
    def getIdentity(cls,):
        if cls.__config['identity'] is None:
            raise Exception('call Configure.setIdentity(\'yourDockerId\') first before using docker')
        return cls.__config['identity']
    
    @classmethod
    def setGenome(cls,val):
        if cls.__config['refdir'] is None:
            raise Exception('refdir should be configure first, call Configure.setRefDir for configuration')
        ## bowtie2 index
        cls.__config['bt2Idx'] = os.path.join(cls.getRefDir(),val)
        suffix = ['.1.bt2','.2.bt2','.3.bt2','.4.bt2','.rev.1.bt2','.rev.2.bt2']
        cls.__config['bt2IdxFiles'] = [ cls.__config['bt2Idx'] + s for s in suffix ]
        
        cls.__config['genome'] = val
        
        
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
        if not folderPath.startswith(cls.getDockerPath()):
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
        cls.startDocker('V1')                      
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
   
    
    @classmethod
    def _getContainerName(cls,version):
        identity = Configure.getIdentity()
        if identity is None:
            identity = ''
        else:
            identity = '_' + identity
        return 'hca_'+ version + identity
    
    @classmethod
    def startDocker(cls, versions = None):        
        if versions is None:
            versions = Configure.getDockerVersions()
        elif not isinstance(versions,list):
            versions = [versions]
        for version in versions:
            try:
                print(' '.join([
                        'docker',
                        'run',
                        '-d',
                        '-t',
                        '--name=' + cls._getContainerName(version),
                        '--privileged=true',
                        '-v',Configure.getTmpDir() + ':' + Configure.getDockerPath(),
                        Configure.getDockerVersion(version),
                        '/bin/bash',
                        ]))
                subprocess.run([
                        'docker',
                        'run',
                        '-d',
                        '-t',
                        '--name=' + cls._getContainerName(version),
                        '--privileged=true',
                        '-v',Configure.getTmpDir() + ':' + Configure.getDockerPath(),
                        Configure.getDockerVersion(version),
                        '/bin/bash',
                        ],stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            except:
                pass
        
    @classmethod
    def stopDocker(cls, versions = None):          
        if versions is None:
            versions = Configure.getDockerVersions()
        elif not isinstance(versions,list):
            versions = [versions]
        for version in versions:
            try:
                subprocess.run([
                        'docker',
                        'rm', 
                        '-f',
                        cls._getContainerName(version),
                        ]) 
            except:
                pass 
    
    @classmethod
    def getDockerCMD(cls, cmdline, version):
        return 'docker exec -it '+ cls._getContainerName(version) +' '+ cmdline
            
    
        

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
        self.__virtual = False
        self.__virtualCount = -1
        self.unsetParams = ''
        self.funGroup = None
        self.__isFinished = False
        self.__multiRun = False
        self.checkParamValid()
        self.__stepID = StepBase.regStepID()
        self.__inputSize = -1
        self.__upstreamSize = 1
        self.tmpdirStack = []
        return 0
    
    def getStepID(self,):
        return self.__stepID
    
   
    def getStepFolderName(self,):
        return 'step_' + str(self.getStepID()).zfill(2) + '_' + self.__class__.__name__
    
    def getStepFolerPath(self,):
        if self.top() is None:
            return Configure.getTmpPath(self.getStepFolderName())
        else:
            return Configure.getTmpDir()
    
    def initIO(self,): 
        if not os.path.exists(self.getStepFolerPath()):
            os.mkdir(self.getStepFolerPath())
        #os.makedirs(os.path.join(Configure.getTmpDir(),self.getStepFolderName(),'.tmp_for_docker',self.getStepFolerPath()))    
        self.push(self.getStepFolerPath())
        #Configure.setTmpDir(self.getStepFolerPath())
        self.impInitIO()
        self.pop()
        #Configure.setTmpDir(tmpdir)
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
    
    def linkRecursive(self,originPath,desPath):
        #print(originPath)
        #print('To')
        #print(desPath)
        if os.path.isfile(originPath):
            subprocess.run(['rm','-rf',desPath])
            subprocess.run(['ln',originPath,desPath])
            return
        elif os.path.isdir(originPath):
            os.makedirs(desPath,exist_ok=True)
            files = os.listdir(originPath)
            for f in files:
                self.linkRecursive(os.path.join(originPath,f),os.path.join(desPath,f))
        else:
            print(['originPath:',originPath,'is neither dir nor file'])
            #raise Exception('originPath:',originPath,'is neither dir nor file')
            
                    
    
    def getIOtype(self,value):
        if value is None:
            return None
        value = os.path.abspath(value)
        for key in self.inputs.keys():
            for path in self.convertToList(self.inputs[key]):
                apath = os.path.abspath(path)
                if apath.startswith(value):
                    if os.path.isdir(value): 
                        return 'inputDir'
                    elif os.path.isfile(value):
                        return 'inputFile'
                    elif os.path.isdir(os.path.dirname(value)):
                        return 'inputPrefix'
        for key in self.outputs.keys():
            for path in self.convertToList(self.outputs[key]):
                apath = os.path.abspath(path)
                if apath.startswith(value):                    
                    if value == apath:
                        return 'outputFile'
        for key in self.outputs.keys():
            for path in self.convertToList(self.outputs[key]):
                apath = os.path.abspath(path)
                if apath.startswith(value):                    
                    if apath[0:(len(value)+1)] == os.path.sep:
                        return 'outputDir'
                    else:                        
                        return 'outputPrefix'
        return None
    
    def linkVirtualPaths(self,origin):
        if not os.path.exists(origin):
            raise Exception('origin path:',origin,'does not exist')
        virDir = os.path.join(
                Configure.getTmpDir(),
                self.getStepFolderName(), 
                '.tmp_for_docker',
                 os.path.dirname(origin)[1:])
        virPath = os.path.join(
                Configure.getTmpDir(),  
                self.getStepFolderName(),
                '.tmp_for_docker',
                 origin[1:])
        print(virDir)
        os.makedirs(virDir,exist_ok=True)
        print(['ln','-f',origin,virPath])
        subprocess.run(['ln','-f',origin,virPath])
    
    def linkRealPaths(self,des):
        virPath = os.path.join(
                Configure.getTmpDir(),  
                self.getStepFolderName(),
                '.tmp_for_docker',
                 des[1:])
        os.makedirs(os.path.dirname(des),exist_ok=True)
        #print(['ln','-f',virPath,des])
        if os.path.exists(virPath):
            subprocess.run(['ln','-f',virPath,des])
        else:
            raise Exception(virPath,des,'does not create succesfull')
        
                
    
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
            
            tmpFolder = Configure.getTmpPath(os.path.join(self.getStepFolderName(),'.tmp_for_docker'))
            os.makedirs(tmpFolder,exist_ok=True)
            print(os.path.join(tmpFolder,'*'))
            subprocess.run(['rm','-rf',os.path.join(tmpFolder,'*')])
            
            if Configure.isDocker(): 
                for key in self.getInputs():
                    files = self.convertToList(self.inputs[key])
                    if files[0] is None:
                        continue
                    for afile in files:
                        self.linkVirtualPaths(afile)
                self.__virtual = True        
                self.push(os.path.join(Configure.getDockerPath(),self.getStepFolderName()))            
            else:
                self.push(self.getStepFolerPath())
           
            
            if self.__multiRun :
                self._multiRun()
                if Configure.isDocker():
                    self.callCmdline('V1', 'chmod 777 -R ' + Configure.getTmpDir(),shell = True, stdoutToLog = False)
                else:
                    pass#self.callCmdline('V1', 'chmod 777 -R ' + self.top(), shell = True, stdoutToLog = False)
            else:
                if self.__inputSize == -1:
                    raise Exception('call self._setInputSize(your sample size) in impInitIO')
                for i in range(self.__inputSize):
                    self._singleRun(i)
                    if Configure.isDocker():
                        self.callCmdline('V1', 'chmod 777 -R ' + Configure.getTmpDir(),shell = True, stdoutToLog = False)
                    else:
                        pass#self.callCmdline('V1', 'chmod 777 -R ' + self.top(), shell = True, stdoutToLog = False)                  
            
            self.pop()        
            if Configure.isDocker():                 
                self.__virtual = False
                for key in self.getOutputs():
                    files = self.convertToList(self.outputs[key])
                    print('1111111111111')
                    print(files)
                    if files[0] is None:
                        continue
                    for afile in files:
                        print(10000000000)
                        print(afile)
                        self.linkRealPaths(afile)


            subprocess.run(['rm','-rf',os.path.join(tmpFolder,'*')])
            
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
                

    
    def checkResult(self,):
        return True
    
    def checkParamValid(self,):
        return True
    
    def getInputs(self,):
        return list(self.inputs.keys())
           
    def getInput(self, inputName):
        if self.__virtual:
            if isinstance(self.inputs[inputName],list):
                #print([os.path.join(Configure.getDockerPath(),self.getStepFolderName(),'.tmp_for_docker',s[1:]) for s in self.inputs[inputName]])
                return [os.path.join(Configure.getDockerPath(),self.getStepFolderName(),'.tmp_for_docker',s[1:]) for s in self.inputs[inputName]]
            else:
                #print(os.path.join(Configure.getDockerPath(),self.getStepFolderName(),'.tmp_for_docker',self.inputs[inputName][1:]))
                return os.path.join(Configure.getDockerPath(),self.getStepFolderName(),'.tmp_for_docker',self.inputs[inputName][1:])
        else:
            return self.inputs[inputName]
    
    def absolutePath(self,pathOrPathList):
        if pathOrPathList is None:
            return None
        elif isinstance(pathOrPathList,list):
            return [os.path.abspath(s) for s in pathOrPathList]
        else:
            return os.path.abspath(pathOrPathList)

    def setInput(self, inputName, inputValue):
        self.inputs[inputName] = self.absolutePath(inputValue)
        if self.inputs[inputName] is not None:
            if isinstance(self.inputs[inputName],list):
                for s in self.inputs[inputName]:
                    os.makedirs(os.path.dirname(os.path.join(Configure.getTmpDir(),'.tmp_for_docker',s[1:])), exist_ok=True)                
                    
            else:
                os.makedirs(os.path.dirname(os.path.join(Configure.getTmpDir(),'.tmp_for_docker',self.inputs[inputName][1:])), exist_ok=True)
                
    
    def getOutputs(self,):
        return list(self.outputs.keys())
    
    def getOutput(self,outputName):
        if self.__virtual:
            if isinstance(self.outputs[outputName],list):
                return [os.path.join(Configure.getDockerPath(),self.getStepFolderName(),'.tmp_for_docker',s[1:]) for s in self.outputs[outputName]]
            else:
                return os.path.join(Configure.getDockerPath(),self.getStepFolderName(),'.tmp_for_docker',self.outputs[outputName][1:])
        else:
            return self.outputs[outputName]
    

    
    def setOutput(self, outputName, outputValue):
        self.outputs[outputName] = self.absolutePath(outputValue)
        if self.outputs[outputName] is not None:
            if isinstance(self.outputs[outputName],list):
                for s in self.outputs[outputName]:
                    os.makedirs(os.path.dirname(os.path.join(Configure.getTmpDir(),
                                                             '.tmp_for_docker',
                                                             s[1:])), exist_ok=True) 
            else:
                os.makedirs(os.path.dirname(os.path.join(Configure.getTmpDir(),'.tmp_for_docker',self.outputs[outputName][1:])), exist_ok=True)
            
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
        if paramName in self.paramsIO.keys():
            if self.__virtual :
                if isinstance(self.paramsIO[paramName],list):
                    return [os.path.join(Configure.getDockerPath(),self.getStepFolderName(),'.tmp_for_docker',s[1:]) for s in self.paramsIO[paramName]]
                else:
                    return os.path.join(Configure.getDockerPath(),self.getStepFolderName(),'.tmp_for_docker',self.paramsIO[paramName][1:])
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
        
    def _setUpstreamSize(self,size):
        self.__upstreamSize = size
    def _setVirtual(self, virtual = True):
        self.__virtual = virtual
    def push(self,newTmp):        
        self.tmpdirStack.append(Configure.getTmpDir())
        Configure.setTmpDir(newTmp)
        return Configure.getTmpDir()
    def pop(self,):
        tpdir=self.tmpdirStack.pop()
        Configure.setTmpDir(tpdir)
        return tpdir
         
    def top(self,):
        if len(self.tmpdirStack) == 0:
            #print('return here')
            return None
        else:
            # print('return here1')
            return self.tmpdirStack[-1]
    
    
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
    
    def getFileList(self, pathPrefix, filePath):
        """
        For developer:
        """
        fileList = self.getListInFile(filePath)
        for i in range(len(fileList)):
            fileList[i] =  os.path.join(pathPrefix,fileList[i])
        return fileList
    
    def getListInFile(self, filePath):
        """
        For developer:
        """
        listFile = open(filePath)
        fileList = listFile.readlines()
        if len(fileList) == 0:
            raise Exception('file',filePath,'can not be empty')
        return fileList
    
    def setOutputDir1To1(self, outputName, outputDir, outputPrefix, outputSuffix, inputName, sep = '.'):
        """
        For developer:
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
        suffixDot = '.'
        if outputSuffix == '':
            suffixDot =''
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
                    outputList.append(prefixs[i] + sep + outputSuffix)
            else:
                if len(inputList) == 1:
                    outputList.append(outputPrefix + suffixDot + outputSuffix)
                else:
                    for i in range(len(inputList)):
                        outputList.append(outputPrefix + sep + str(i) + suffixDot + outputSuffix)
            if outputDir is None:
                self.setOutput(outputName,Configure.getTmpPath(outputList))
            else:        
                self.setOutput(outputName,[os.path.join(outputDir,s) for s in outputList])
                
    def setOutputDir1To1ByFunc(self, outputName, outputDir, func, inputName):
        """
        For developer:
        Use the known input file paths to generate the output file paths.
        outputName(str): the key name of the file path list to be generate
        outputDir(str or None): the folder of the output file stored, None or string is OK.
                   if None, its value will set globally
        func: an function of func(fileName),it shoul implement: the fileName are provided and a generated output file name should be return)
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
            fileBasename = [os.path.basename(s) for s in inputList]
            for i in range(len(inputList)):
                outputFile = func(fileBasename[i])
                if not isinstance(outputFile,str):
                    raise Exception('file name generated by func is not a string')
                outputList.append(outputFile)            
            if outputDir is None:
                self.setOutput(outputName,Configure.getTmpPath(outputList))
            else:        
                self.setOutput(outputName,[os.path.join(outputDir,s) for s in outputList])
                
    def setOutputDirNTo1(self, outputName, outputFilePath, fileNameForUnsetOutPath, inputName):
        """
        For developer:
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
                self.setOutput(outputName,outputFilePath)
            
    def callCmdline(self,dockerVersion,cmdline,shell = False, stdoutToLog = True, otherPrefix = None):
        """
        For developer:
            call the command line and write log 
        """
        
        if not shell:
            cmdline = ' '.join(cmdline)
            
        if Configure.isDocker() and dockerVersion is not None:
            cmdline = 'bash -c \"' + cmdline +'\"'            
            if self.top() is not None:
                
                cur = Configure.getTmpDir()
                self.pop()                
                Schedule.startDocker(versions=dockerVersion)
                self.push(cur)
            else:                
                Schedule.startDocker(versions=dockerVersion) 
                
            cmdline = Schedule.getDockerCMD(cmdline,dockerVersion)
        
        if otherPrefix is not None:
            cmdline = cmdline + ' ' + otherPrefix
        
        
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
        """
        For developer:
            get command line parameters for boolean parameter
        """
        return cmdName if self.getParam(paramName) else ''
    
    def setInputFileInDir(self, inputName, inputDir, inputFileName):
        """
        For developer:
        when input parameter is a list or single string,
        use this parameter to generate all of the file paths under the path or
        just use the single string path directly
        inputName: the key name of inputs
        inputValue: string of directory or file path, or list of file paths  
        """
        if inputDir is None:
            self.setInput(inputName, None)
        else:
            if inputFileName is None:
                raise Exception('inputFileName can not be None')
            else:
                self.setInput(inputName,os.path.join(inputDir,inputFileName))
                    
    def setInputDirOrFile(self, inputName, inputValue, checkSameSuffix=True):
        """
        For developer:
        when input parameter is a list or single string,
        use this parameter to generate all of the file paths under the path or
        just use the single string path directly
        inputName: the key name of inputs
        inputValue: string of directory or file path, or list of file paths  
        """
        
        if inputValue is None:
            self.setInput(inputName, None)
        else:
            if isinstance(inputValue,list):
                self.setInput(inputName, inputValue)
            else:
                if os.path.isdir(inputValue):
                    filelist = os.listdir(inputValue)
                    filelist.sort()
                    suffix = [s.split('.')[-1] for s in filelist]
                    if checkSameSuffix and len(set(suffix)) != 1:
                        raise Exception('the suffix of files under path:',inputName,'is not the same, check the file format under the directory')
                    self.setInput(inputName, [os.path.join(inputValue,s) for s in filelist ])
                else:
                    self.setInput(inputName, inputValue)
                    if os.path.exists('/home/zhengwei/zuoye/step_00_AdapterRemoval/step_00_AdapterRemoval/'):
                        raise Exception('7777777777777777777777777777777777777777')
        
                    
    def setInputDirRecursion(self,inputValue):
        self.setInputDirRecursionFunc(inputValue,'')
        
    def setInputDirRecursionFunc(self,inputValue, paFolder):
        if inputValue is None:
            return 
        if not os.path.isdir(inputValue):
            raise Exception('inputValue is not a directory:',inputValue)
        fileOrDir = os.listdir(inputValue)
        fileOrDir.sort()
        
        fileparFoder = [os.path.join(paFolder,s) for s in fileOrDir]
        fileOrDir = [os.path.join(inputValue,s) for s in fileOrDir]        
        for i in range(len(fileOrDir)):
            if os.path.isdir(fileOrDir[i]):
                self.setInputDirRecursionFunc(fileOrDir[i],fileparFoder[i])
            else:
                self.setInput(fileparFoder[i],fileOrDir[i])
      
    
    def getInputList(self, inputName):
        """
        For developer:
        the wraper of getInput(). if the input is the single string,
        it will wrap it will list. Otherwise it will return the input as it is.
        """
        return self.convertToList(self.getInput(inputName))         

    def getInputDir(self, inputName):
        filePath = self.getInputList(inputName)
        return os.path.dirname(filePath[0])
    
    def getInputPrefix(self, inputName, prefix):       
        return self.getInputDir(inputName) + prefix
        
    def getOutputList(self, outputName):
        """
        For developer:
        the wraper of getOutput(). if the output is the single string,
        it will wrap it will list. Otherwise it will return the output as it is.
        """
        return self.convertToList(self.getOutput(outputName))

    def getOutputDir(self, outputName):
        filePath = self.getOutputList(outputName)
        return os.path.dirname(filePath[0])
    
    def getOutputPrefix(self, outputName, prefix):
        return self.getOutputDir(outputName) + prefix
    
    
    def getInputs(self,):
        """
        For developer and user
        get Valid keys for input file paths 
        """
        return super(Step,self).getInputs()
           
    def getInput(self, inputName):
        """
        For developer and user
        get input file paths of inputName
        """
        return super(Step,self).getInput(inputName)

    def setInput(self, inputName, inputValue):
        """
        For developer and user
        set input file paths of inputName with inputValue
        """
        super(Step,self).setInput(inputName,inputValue)
    
    def getOutputs(self,):
        """
        For developer and user
        get Valid keys for output file paths 
        """
        return super(Step,self).getOutputs()
    
    def getOutput(self,outputName):
        """
        For developer and user
        get output file paths of outputName
        """
        return super(Step,self).getOutput(outputName)
    
    def setOutput(self, outputName, outputValue):
        """
        For developer
        set output file paths of outputName with inputValue
        """
        super(Step,self).setOutput(outputName,outputValue)
    
    def getUnsetParams(self,):
        return super(Step,self).getUnsetParams()
    
    def getParams(self,):
        """
        For developer and user
        get Valid keys for parameters except for input and output parameters 
        """
        return super(Step,self).getParams()
    
    def getParam(self, paramName):
        """
        For developer and user
        get parameter of paramName. parameters except for input and output parameters
        """
        return super(Step,self).getParam(paramName)
    
    def setParam(self, paramName, paramValue):
        """
        For developer
        Set parameters (except for input or output parameters) in __init__()
        """
        super(Step,self).setParam(paramName, paramValue)
        
    def getParamIOs(self,):
        """
        For developer and user
        Get parameters (except for input or output parameters) setted in __init__()
        """
        return super(Step,self).getParamIOs()
    
    def getParamIO(self, paramName):
        """
        For developer and user
        Get input or output parameters set in __init__() or call()
        """
        return super(Step,self).getParamIO(paramName)
    
    def setParamIO(self, paramName, paramValue, setTmp = False):
        """
        For developer
        Set input  or output parameters at __init__() and 
        Set input parameters from upstream at call()
        """

        return super(Step,self).setParamIO(paramName, paramValue)
        
    def initParam(self,**kwargs):
        """
        For User
        reinit the all or part of parameters by calling like this:
        object.initParam(paramName1,'paramValue1',paramName1,'paramValue1')
        """
        keyIO = self.getParamIOs()
        keyParam = self.getParams()
        isReinitIO = False
        for key in list(kwargs.keys()):
            if key not in keyIO and key not in keyParam:
                raise Exception('parameter name:',key,'is not valid')
            
        for key in kwargs.keys():
            if key in keyIO:
                self.setParamIO(kwargs[key])
                isReinitIO = True
            elif key in keyParam:
                self.setParam(kwargs[key])
        
        if isReinitIO:
            self.initIO()
            
    def convertToRealPath(self, virtualPath): 
        if Configure.isDocker():
            if isinstance(virtualPath,list):
                return [os.path.join(self.top(),s[len(Configure.getDockerPath())+1:]) for s in virtualPath]
            else:   
                return os.path.join(self.top(),virtualPath[len(Configure.getDockerPath())+1:])
        else:
            return virtualPath
        
    def convertToRealRealPath(self, virtualPath): 
        if Configure.isDocker():
            if isinstance(virtualPath,list):
                return [os.path.join(self.top(),s.split('.tmp_for_docker')[1][1:]) for s in virtualPath]
            else:   
                return os.path.join(self.top(),virtualPath.split('.tmp_for_docker')[1][1:])
        else:
            return virtualPath
    
    def getMarkdown(self,lang='EN'):        
        if self.checkFinish():
            self.push(self.getStepFolerPath())
            if lang == 'EN':
                return self.getMarkdownEN()
            elif lang == 'CN':
                return self.getMarkdownCN()
            else:
                raise Exception('language',lang,'is not support yet!')
            self.pop()
        else:
            raise Exception(self.getStepFolderName(),'is not finished')
    
    def getMarkdownEN():
        raise Exception('getMarkdownEN must be overwrote')
        
    def getMarkdownCN():
        raise Exception('getMarkdownCN must be overwrote')
        
        