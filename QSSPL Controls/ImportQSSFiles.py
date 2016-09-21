
import numpy as np
import os
from shutil import copyfile

class Load_QSSPL_File_LabView():
    Directory  = ''
    RawDataFile = ''

    def Load_RawData_File_LabView(self):
        return np.genfromtxt(self.Directory+self.RawDataFile,names= ('Time','PC','Gen','PL'))

    def Load_InfData_File_LabView(self):
        InfFile = self.RawDataFile[:-13]+'.inf'

        '''info from inf file '''

        Cycles,dump,Frequency,LED_Voltage,dump,dump,dump,dump,DataPoints,dump = np.genfromtxt(self.Directory+InfFile,skip_header=20,skip_footer=22,delimiter=':',usecols=(1),autostrip=True,unpack=True)
        Waveform,LED_intensity     = np.genfromtxt(self.Directory+InfFile,skip_header=31,skip_footer=20,delimiter=':',usecols=(1),dtype=None,autostrip=True,unpack=True) 


        l = np.genfromtxt(self.Directory+InfFile,skip_header=36,delimiter=':',usecols=(1))



        Doping = l[9]
        Ai = l[6]
        Fs = l[7]
        Thickness = l[12]
        Quad = l[12]
        Lin = l[12]
        Const = 0

        Binning = int(l[2])
        Reflection = (1-l[16])*100

        dic = locals()

        del dic['self']
        del dic['l']
        del dic['dump']

        return dic

    def Load_ProcessedData_File_LabView(self):  
        DataFile = self.DataFile[:-13]+'.dat'
        return np.genfromtxt(self.Directory+DataFile,usecols=(0,1,8,9),unpack=True,delimiter='\t',names=('Deltan_PC','Tau_PC','Deltan_PL','Tau_PL'))


    def WriteTo_Inf_File_LabView(self,Dictionary):



        InfFile = self.RawDataFile[:-13]+'.inf'

        if (os.path.isfile(self.Directory+InfFile+".Backup")==False):
            copyfile(self.Directory+InfFile,self.Directory+InfFile+".Backup")
            print 'Backuped original .inf  file as .inf.backup'

        ####
        # Creating the .inf file, this can be done more easily with list(f), but i'm not using it right now.
        ###
        f = open(self.Directory+InfFile, 'r')

        #print list(f)
        #print list(f).shape

        s=''
        for i in range(38):
            s= s+ f.readline()
        s= s+f.readline()[:26]+'{0:.0f}'.format(Dictionary['Binning'])+'\n'
        for i in range(3):
            s= s+ f.readline()
        s= s+f.readline()[:5]+'{0:.3e}'.format(Dictionary['Ai'])+'\n'
        s= s+f.readline()[:11]+'{0:.3e}'.format(Dictionary['Fs'])+'\n'
        s= s+f.readline()
        s= s+f.readline()[:23]+'{0:.3e}'.format(Dictionary['Doping'])+'\n'
        s= s+f.readline()
        s= s+f.readline()
        s= s+f.readline()[:12]+'{0:.4f}'.format(Dictionary['Thickness'])+'\n'
        s= s+f.readline()[:24]+'{0:.4e}'.format(Dictionary['Quad'])+'\n'
        s= s+f.readline()[:21]+'{0:.4e}'.format(Dictionary['Lin'])+'\n'
        s= s+f.readline()
        s= s+f.readline()[:37]+'{0:.6f}'.format(1-Dictionary['Reflection']/100)+'\n'

        for i in range(6):
            s= s+ f.readline()

        f.close()
        f= open(self.Directory+InfFile, 'w')
        f.write(s)









class Load_QSSPL_File_Python():
    Directory  = ''
    RawDataFile = ''
    def Load_RawData_File_Python(self):    
        data = np.genfromtxt(self.Directory+self.RawDataFile,unpack=True,names= True,delimiter='\t')
        s=np.array([])
        dic = {'Time_s':'Time','Generation_V':'Gen','PL_V':'PL','PC_V':'PC'}
        # print np.array(data.dtype.names)
        for i in np.array(data.dtype.names):
            # print i,dic[i]
            s= np.append(s,dic[i])


        # print s

        data.dtype.names = s
        # ('Time','Gen','PL','PC')
        return data

    def num(self,s):
        try:
            return float(s)
        except ValueError:
            return s

    def Load_InfData_File_Python(self):
        # print 'Still under construction'

        InfFile = self.RawDataFile[:-13]+'.inf'
        
        



        #These are adjustment Values
        Doping = 1
        Thickness = 1
        Binning = 1
        Reflection = 0.0
        Fs = 1
        Ai = 1
        Quad = 0.0004338
        Lin = 0.03611



        CropStart = 0
        CropEnd = 100


        List =  locals()
        
        del List['InfFile']
        del List['self']


        with open(self.Directory+str(InfFile), 'r') as f:
            s=f.read()

        for i in s.split('\n')[2:-1]:
            # print i.split(':\t')[1]
            List[i.split(':\t')[0].strip()] = self.num(i.split(':\t')[1])
        # print List
        

        return List
        

    def Load_ProcessedData_File_Python(self):
        print 'Still under construction'

 
        return zeros(4,4)


    def WriteTo_Inf_File_Python(self,Dictionary):

        InfFile = self.RawDataFile[:-13]+'.inf'

        if (os.path.isfile(self.Directory+InfFile+".Backup")==False):
            copyfile(self.Directory+InfFile,self.Directory+InfFile+".Backup")
            print 'Backuped original .inf  file as .inf.backup'


        s='MJ system\r\nList of vaiables:\r\n'
        for i in Dictionary:
            # if i != and i !='self':
                s += '{0}:\t{1}\r\n'.format(i,Dictionary[i])

        
        with open(self.Directory+InfFile,'w') as text_file:
               text_file.write(s)











class LoadData(Load_QSSPL_File_Python,Load_QSSPL_File_LabView):

    Directory  = ''
    RawDataFile = ''
    File_Type = ''

    def Determine_File_Type(self):
        if '.Raw Data.dat' in self.RawDataFile:
            self.File_Type = 'Python'
        elif '_Raw Data.dat' in self.RawDataFile:
            self.File_Type = 'LabView'
        else:
            print 'Unknown file extension'

    def Load_RawData_File(self):
        self.Determine_File_Type()
        return  getattr(self,'Load_RawData_File_'+self.File_Type)()

    def Load_InfData_File(self):
        self.Determine_File_Type()
        return  getattr(self,'Load_InfData_File_'+self.File_Type)()
        
    def Load_ProcessedData_File(self):
        self.Determine_File_Type()
        return  getattr(self,'Load_ProcessedData_File_'+self.File_Type)()

    def WriteTo_Inf_File(self,Dict):


        return  getattr(self,'WriteTo_Inf_File_'+self.File_Type)(Dict)     

if __name__ == "__main__":

    Load = LoadData()



    Load.Directory = r'C:\Users\mattias\Desktop'
    Load.RawDataFile = r'\a.Raw Data.dat'

    Data= Load.Load_RawData_File()
    # print len(Data.dtype.names)
    # print np.mean(Data,axis=0)
    # print Load.WriteTo_Inf_File(Data)
    Load.Load_InfData_File()
    # print Load.Load_RawData_File()['Time'].shape