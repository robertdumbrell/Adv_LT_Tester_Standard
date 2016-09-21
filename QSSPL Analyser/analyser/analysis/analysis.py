

from numpy import *
import matplotlib.pylab as plt

from models.ConstantsClass import *

from models.NumericalDifferentiation_windows import Finite_Difference, Regularisation

from utils.importexport import LoadData
import scipy.constants as C

import caculate_quantities as CQ


def find_nearest(array, value):
    idx = (abs(array - value)).argmin()
    return idx


class Data(Constants):

    Derivitive = 'Regularised'
    Analysis = 'Generalised'
    Type = 'p'
    CropStart = 0
    CropEnd = 100
    BackgroundSubtraction = 0.95
    Used = False
    Temp = -1

    def __init__(self):
        self.LD = LoadData()

    def _update_ni(self, model_handeller):
        if self.Temp != self.Wafer['Temp']:
            self.ni = model_handeller.update['ni'](temp=self.Wafer['Temp'])
            self.Temp = self.Wafer['Temp']
            self.Vt = C.k * self.Temp / C.e
        pass

    def BackgrounConcentration(self):
        if (self.Type == 'p'):
            self.nh0 = self.Wafer['Doping']
            self.ne0 = self.ni**2 / self.Wafer['Doping']

        elif(self.Type == 'n'):
            self.ne0 = self.Wafer['Doping']
            self.nh0 = self.ni**2 / self.Wafer['Doping']
        else:
            self.ne0 = 1e20
            self.nh0 = 1e20

        pass

    def ProvideRawDataFile(self, Directory, RawDataFile):
        self.Directory = Directory + '/'
        self.RawDataFile = RawDataFile

        self.Used = True

        self.Load_Measurements()

    def UpdateInfData(self):
        self.LD.WriteTo_Inf_File(self.Wafer)

    def dndt(self, Deltan):

        if (self.Derivitive == 'Regularised'):
            dn_dt = Regularisation().FirstDerivative(
                self.Data['Time'], Deltan, 1e-20)

        elif (self.Derivitive == 'Finite Difference'):
            dn_dt = Finite_Difference().FourPointCentral(
                self.Data['Time'], Deltan)

        else:
            print('You fucked up.... again')

        return dn_dt

    def Load_Measurements(self):

        self.LD.Directory = self.Directory
        self.LD.RawDataFile = self.RawDataFile

        self.RawData = self.LD.Load_RawData_File()
        self.RawData2 = copy(self.RawData)
        self.Wafer = self.LD.Load_InfData_File()
        print(self.Wafer)

    def ChoosingDefultCropValues(self):

        # if no values are provided, go forth a crop
        if self.Wafer['CropStart'] == None and self.Wafer['CropStart'] == None:

            Waveform = self.Wafer['Waveform']

            if Waveform == 'Triangle':
                self.Wafer['CropStart'], self.Wafer['CropEnd'] = 35, 55
            elif Waveform == 'Square':
                self.Wafer['CropStart'], self.Wafer['CropEnd'] = 13, 50
            elif Waveform == 'Sawtooth':
                self.Wafer['CropStart'], self.Wafer['CropEnd'] = 12, 79
            else:
                self.Wafer['CropStart'], self.Wafer['CropEnd'] = 5, 95

    def iVoc(self):
        PC_ivoc = CQ.iVoc_from_carriers(
            self.ne0, self.nh0, self.DeltaN_PC, self.Wafer['Temp'], self.ni)
        PL_ivoc = CQ.iVoc_from_carriers(
            self.ne0, self.nh0, self.DeltaN_PL, self.Wafer['Temp'], self.ni)
        return PC_ivoc, PL_ivoc

    def CalculateLifetime(self, BackGroundShow=False, model_handeller=None):

        # make sure the ni is updated
        self._update_ni(model_handeller)

        # determine the background concentration of carriers
        self.BackgrounConcentration()

        # Background correction stuff
        BackgroundIndex = int(self.RawData['Time'].shape[0] * .05)
        self.Data = copy(self.RawData)

        for i in ['PL', 'Gen']:
            self.Data[i] -= average(self.Data[i][:BackgroundIndex])

        self.Data['PC'] = self.Wafer['Quad'] * \
            (self.Data['PC']) * (self.Data['PC']) + \
            self.Wafer['Lin'] * (self.Data['PC']) + self.Wafer['Const']

        '''for now just use the 95% limit to cal background'''
        """background subtraction"""
        self.DarkConductance = average(self.Data['PC'][:BackgroundIndex])
        for i in ['PC']:

            self.Data[i] -= average(self.Data[i][:BackgroundIndex])

        self.Cropping_Percentage()

        # """For checking background subtraction"""
        # if BackGroundShow == True:
        #     fig = plt.figure('BackGround Check')
        #     plt.title('BackGround Check')
        #     for i in ['PC', 'PL', 'Gen']:
        #         plt.plot(self.Data['Time'], self.Data[i], label=i)
        #     plt.xlim(0, max(self.Raw_Time))
        #     plt.legend(loc=0)
        #     plt.ylim(-1, 11)
        #     plt.show()

        self.Data = self.Binning_Named(self.Data, self.Wafer['Binning'])

        model_handeller._update_update()
        self.DeltaN_PC = CQ.nxc_from_photoconductance(
            self.Data['PC'],
            self.Wafer['Thickness'],
            self.Wafer['Temp'],
            self.ne0,
            self.nh0,
            model_handeller)

        if self.Wafer['Type'] == 'n':
            dopant = 'phosphorous'
        elif self.Wafer['Type'] == 'p':
            dopant = 'boron'

        self.DeltaN_PL = CQ.nxc_from_photoluminescence(
            self.Data['PL'],
            self.Wafer['Ai'],
            dopant,
            abs(self.ne0 - self.nh0),
            self.Wafer['Temp'],
            model_handeller)

        self.Tau_PC = self.DeltaN_PC / self.Generation('PC')
        self.Tau_PL = self.DeltaN_PL / self.Generation('PL')

    def Generation(self, PCorPL, suns=False):
        try:
            Gen = getattr(
                self, 'Generation_' + self.Analysis.replace(' ', '_'))
        except:
            print('Choice of generation doesn\'t exist: You fucked up')

        if suns == True:
            scale = self.Wafer['Thickness'] / 2.5e17
        else:
            scale = 1.

        return Gen(PCorPL) * scale

    def Generation_Steady_State(self, PCorPL):
        Trans = (1 - self.Wafer['Reflection'] / 100.)
        return self.Data['Gen'] * self.Wafer['Fs'] * Trans / self.Wafer['Thickness']

    def Generation_Generalised(self, PCorPL):
        Trans = (1 - self.Wafer['Reflection'] / 100.)
        if PCorPL == 'PC':
            return self.Data['Gen'] * self.Wafer['Fs'] * Trans / self.Wafer['Thickness'] - self.dndt(self.DeltaN_PC)
        elif PCorPL == 'PL':
            return self.Data['Gen'] * self.Wafer['Fs'] * Trans / self.Wafer['Thickness'] - self.dndt(self.DeltaN_PL)
        else:
            print('You fucked up the Generation')

    def Generation_Transient(self, PCorPL):
        if PCorPL == 'PC':
            return -self.dndt(self.DeltaN_PC)
        elif PCorPL == 'PL':
            return -self.dndt(self.DeltaN_PL)
        else:
            print('You fucked up the Generation')

    def Local_IdealityFactor(self):
        # Generation scale doesn't matter so Generation is used
        iVocPC, iVocPL = self.iVoc()

        if (self.Derivitive == 'Regularised'):
            return  (self.Generation('PC')) / (self.Vt * Regularisation().FirstDerivative       (self.Data['Time'], self.Generation('PC'), 1e-20) / Regularisation().FirstDerivative  (self.Data['Time'], iVocPC, 1e-20)),\
                    (self.Generation('PL')) / (self.Vt * Regularisation().FirstDerivative(self.Data['Time'], self.Generation(
                        'PL'), 1e-20) / Regularisation().FirstDerivative(self.Data['Time'], iVocPL, 1e-20))

        elif (self.Derivitive == 'Finite Difference'):
            return  (self.Generation('PC')) / (self.Vt * Finite_Difference().FourPointCentral   (self.Data['Time'], self.Generation('PC'))     / Finite_Difference().FourPointCentral(self.Data['Time'], iVocPC)),\
                    (self.Generation('PL')) / (self.Vt * Finite_Difference().FourPointCentral(self.Data[
                        'Time'], self.Generation('PL')) / Finite_Difference().FourPointCentral(self.Data['Time'], iVocPL))

    def Cropping_negitives(self):
        # this just uses that when points are negitive they should no longer be
        # used
        maxindex = argmax(self.SS_Generation)

        IndexOfNegitives = where(self.SS_Generation < 0)[0]

        firstnegtive = where(IndexOfNegitives < maxindex)[0][-1]
        lastnegtive = where(IndexOfNegitives > maxindex)[0][0]
        self.index = arange(
            IndexOfNegitives[firstnegtive], IndexOfNegitives[lastnegtive], 1)
        # plot(self.Time,self.DeltaN_PL)
        # show()

        # IndexOfNegitives,maxindex,self.DeltaN_PL[IndexOfNegitives[firstnegtive]]

        #self.DeltaN_PL = self.DeltaN_PL [ self.index]
        #self.DeltaN_PC = self.DeltaN_PC [ self.index]
        self.SS_Generation = self.SS_Generation[self.index]
        self.Time = self.Time[self.index]

        self.RawPCDataEdited = self.RawPCDataEdited[self.index]
        self.Raw_PLEdited = self.Raw_PLEdited[self.index]

    def Cropping_Percentage(self):
        # this just uses that when points are negitive they should no longer be
        # used
        maxindex = amax(self.Data.shape)

        self.index = arange(int(self.Wafer[
                            'CropStart'] / 100 * maxindex), int(self.Wafer['CropEnd'] / 100 * maxindex), 1)

        self.Data = self.Data[self.index]

    def EQE(self):

        return self.RawPCDataEdited / self.Generation('PC'), self.Raw_PLEdited / self.Generation('PL') * self.Ai

    def IQE_SingleIntensity(self):

        idx = find_nearest(self.Generation('PL'), 2.5e17 / self.Thickness / 10)
        return self.RawPCDataEdited[idx] / self.Generation('PC')[idx], self.Raw_PLEdited[idx] / self.Generation('PL')[idx] * self.Ai
