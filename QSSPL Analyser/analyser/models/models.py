

from semiconductor.electrical.mobility import Mobility
from semiconductor.electrical.ionisation import Ionisation as Ion
from semiconductor.material.ni import IntrinsicCarrierDensity as NI
from semiconductor.recombination.intrinsic import Radiative


class models_handeller():
    '''

    '''

    available_models = {}
    material = 'Si'

    def __init__(self):
        self._get_available_models()
        self.selected_model = {'ni': 'Couderc_2014',
                               'mobility': 'klaassen_1992',
                               'ionisation': 'Altermatt2006_table1',
                               'B': 'Altermatt2005'
                               }
        self._update_update()

    def access(self, model, material, author, **kwargs):

        return self.use_models[model](
            material=material, author=author
        ).update(kwargs)

    def _get_available_models(self):
        ni = NI().available_models()
        mobility = Mobility().available_models()
        ionisation = Ion().available_models()
        B = Radiative().available_models()

        values = locals()
        del values['self']
        self.available_models = values
        return self.available_models

    def _update_update(self):
        '''
        creates a dictionary that holds
        the required semiconudctor models for easy
        calling
        '''

        self.update = {
            'ni': NI(
                material=self.material,
                author=self.selected_model['ni']
            ).update,
            'mobility': Mobility(
                material=self.material,
                author=self.selected_model['mobility']
            ).mobility_sum,
            'ionisation': Ion(
                material=self.material,
                author=self.selected_model['ionisation']
            ).update_dopant_ionisation,
            'B': Radiative(
                material=self.material,
                author=self.selected_model['B'],
            )._get_B
        }

    def _auto_select_models(self):
        values = self._get_available_models()

        self.selected_model = {}
        for i, k in zip(values.keys(), values.values()):
            self.selected_model[i] = k[0]
