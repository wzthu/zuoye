

from ..core import Step,Configure
from ..steps import SingleCellExperiment
import os
class SC3_DE(Step):
    def __init__(self,
                 sceInput = None,
                 outputpath = None,
                 cluster_num = 0,
                 cmdParam=None,
                 **kwargs):
        """
        SC3_DE is XXX. Need to use this 
        Step as the downstream of SingleCellExperiment.
        >SC3_DE():_init_parameters
            sce: str
            The R workspace saved by SingleCellExperiment Step.
            outputpath: str
            A str indicates the name of appointed folder that saves outputs.You should
            build that folder in advance. The absolute path is also legel. 
            cluster_num: int
            The number of clusters.
            set to 0 will auto estimate the cluster number
            cmdParam: str or list of string
            current unsupported
        >SC3_DE()():_call_parameters
            Avaliabel upstream objects combinations:
            (SingleCellExperiment)
        """
        super(Step, self).__init__(cmdParam,**kwargs)
        
        

        # set all input and output parameters
        self.setParamIO('sceInput',sceInput)
        self.setParamIO('outputpath',outputpath) 
        # call self.initIO()
        self.initIO()
        #set other parameters


        self.setParam('cluster_num',cluster_num)
        #self._setMultiRun()
        
    def impInitIO(self,):
        """
        This function is to initialize 
        all of the input and output files from the io parameters set in __init__() 
        """
        # obtain all input and output parameters        
               
        outputpath = self.getParamIO('outputpath')  
        #set the input file
        if outputpath is None:
            self.setParamIO('outputpath',Configure.getTmpDir()) 
            outputpath = self.getParamIO('outputpath') 

        sceInput = self.getParamIO('sceInput') 
        #set all input files
        self.setInputDirOrFile('sceInput',sceInput)

        # create output file paths and set
        #self.setOutputDir1To1ByFunc('sceOutput',outputpath,func,"matrix_file")
        #self.setOutputDir1To1('sc3OutputFolder',outputpath,None,"_folder","sceInput",sep='')
        def func1(basename):
            return basename+'/Expression.jpg'
        def func2(basename):
            return basename+'/DE_Gene.jpg'
        def func3(basename):
            return basename+'/Gene_Marker.jpg'          
        self.setOutputDir1To1ByFunc('sc3Output_Expression.jpg',outputpath, func1,'sceInput')
        self.setOutputDir1To1ByFunc('sc3Output_DE_Gene.jpg',outputpath, func2,'sceInput')
        self.setOutputDir1To1ByFunc('sc3Output_Gene_Marker.jpg',outputpath, func3,'sceInput')

        # Rscripts
        self.setInputRscript('Rscript','SC3.R')

        if sceInput is not None:
            self._setInputSize(len(self.getInputList('sceInput')))
    def call(self,*args):

        Upstream = args[0]
        if isinstance(Upstream,SingleCellExperiment):
            self.setParamIO('sceInput', Upstream.getOutput('sceOutput'))


    def getMarkdownEN(self,):
        Expression = self.getOutput('sc3Output_Expression.jpg')
        DE_Gene = self.getOutput('sc3Output_DE_Gene.jpg')
        Gene_Marker = self.getOutput('sc3Output_Gene_Marker.jpg')

        Expression_sen = ['***For %s\n***\n![Expression](%s)'%(item.split("/")[-2],item) for item in Expression]
        DE_Gene_sen = ['***For %s\n***\n![DE_Gene](%s)'%(item.split("/")[-2],item)for item in DE_Gene]
        Gene_Marker_sen = ['***For %s\n***\n![Gene_Marker](%s)'%(item.split("/")[-2],item)for item in Gene_Marker]
        Expression_sen = "\n".join(Expression_sen)
        DE_Gene_sen = "\n".join(DE_Gene_sen)
        Gene_Marker_sen = "\n".join(Gene_Marker_sen)
        mdtext = """
## SC3_DE Usage

SC3_DE('/path/to/sce.RData','/path/to/output_dir',cluster_num)  

## SC3 Differential Expression Result  
The SC3 Differential Expression result is shown below:  

### Expression Result  

{Expression} 
### Differential Expression Detecet Result  

{DE_Gene} 

### Gene Marker Detect Result  

{Gene_Marker}     
""".format(Expression=Expression_sen ,
    DE_Gene =DE_Gene_sen,
    Gene_Marker =Gene_Marker_sen,
        )


        return mdtext
            
            
    def _singleRun(self,i):
        # obtain all input and output dir list
        sceInputs = self.getInputList('sceInput')
        cluster_num =  self.getParam('cluster_num')
        sc3Outputjpgs = self.getOutput('sc3Output_Expression.jpg')
        os.path.dirname
        Rscript = self.getInput('Rscript')
        cmdline =['Rscript',
                  Rscript,
                   sceInputs[i],
                   str(cluster_num),
                   os.path.dirname(sc3Outputjpgs[i]),
                   'de'
                   ]
        self.callCmdline('V1', cmdline)