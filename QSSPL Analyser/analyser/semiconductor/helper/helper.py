#!/usr/local/bin/python
# UTF-8

import matplotlib.pylab as plt
import numpy as np
import json
import inspect
try:
    import ConfigParser as configparser
except:
    import configparser


def change_model(Models, author=None):

    author = author or Models.get('default', 'model')

    model = Models.get(author, 'model')

    vals = dict(Models.items(author))

    del vals['model']

    for k, v in vals.items():

        # checks if float or list
        try:
            vals[k] = float(v)
        except:
            try:
                vals[k] = [float(i) for i in v.split(';')]
            except:
                pass

    return vals, model


class HelperFunctions(object):

    _cal_dts = {
        'material': 'Si',
        'temp': 300,
    }

    def __init__(self):
        pass

    @property
    def calculationdetails(self):
        return self._cal_dts

    @calculationdetails.setter
    def calculationdetails(self, kwargs):
        '''
        assignes the inputted values that are requrired,
        befor calling a function to pass it to the downstream
        classes
        '''
        if kwargs:

            items = [i for i in kwargs.keys() if i in self._cal_dts.keys()]

            for item in items:
                if item in self._cal_dts.keys():
                    self._cal_dts[item] = kwargs[item]
                else:
                    print('item:', item)
                    print('not in dictionary:', self._cal_dts.keys())
                    print(self._cal_dts['author'])

    #         self._update_links()
    # def _update_links(self):
    #     pass

    def _int_model(self, fname):
        self.Models = configparser.ConfigParser()
        self.Models.read(fname)

    def change_model(self, author, Models=None):

        Models = Models or self.Models

        self.vals, self.model = change_model(Models, author)

    def plot_all_models(self, update_function, xvalues=None, **kwargs):
        '''
        cycles through all the models and plots the result
        inputs:
            update_function: str
                 the name of the specific function used to update the author
                i.e 'update_Eg'
            **kwargs:
                variables to be passed to the update function.
        '''
        fig, ax = plt.subplots(1)
        for model in self.available_models():

            self.change_model(model)
            result = getattr(self, update_function)(**kwargs)
            if xvalues is None:
                ax.plot(result, label=model)
            else:
                ax.plot(xvalues, result, label=model)
        plt.legend(loc=0)

    def plotting_colours(self, n_colours, fig, ax, repeats=None):
        '''
        a function to help with plotting
        Lets you get a unique range of colours,
        and have repeats of colours
        '''

        colours = []

        # Define the colours for the figures, using CMRmap to get the colours
        for i in np.linspace(.1, .7, n_colours):

            # if don't want any repeat of colours
            if repeats is None:
                colours.append(plt.cm.CMRmap(i))
            else:
                for j in range(repeats):
                    colours.append(plt.cm.CMRmap(i))

        # incase its a 2D axis
        try:
            for axes in ax.flatten():
                axes.set_color_cycle(colours)
        except:
            ax.set_color_cycle(colours)

    def available_models(self, Filter=None, Filter_value=None):
        '''
        Returns a list of all authors for a given parameter this a class of.
        The returned list can be filtered in field "Filter" by the value
        "Filter_value"

        inputs: (all optional)
            Filter: (str)
                The model field that is to be checked
            Filter_value:
                The value to be checked for
        returns:
            list of authors
        '''

        author_list = self.Models.sections()
        author_list.remove('default')

        # remove modles that are not implimented
        for author in list(author_list):

            if 'not_implimented' in dict(self.Models.items(author))['model']:
                author_list.remove(author)

        # does the filtering
        if Filter is not None:
            for c, author in enumerate(author_list):
                if self.Models.get(author, Filter) not in Filter_value:
                    author_list.remove(author)

        # prints no models available
        if not author_list:
            print('No authors for this models available')

        return author_list

    def print_model_notes(self, model=None):
        '''
         prints the notes about the modells
         inputs:
            model: str
                prints notes on model, if not model is seltecte prints
                all model notes
        '''

        if model is None:
            models = self.available_models()
        else:
            models = [model]

        for mdl in models:
            print('{0}:\t'.format(mdl),)
            try:
                # print
                print(dict(self.Models.items(mdl))['notes'])
            except:
                print('No notes')


class Webplotdig_JSONreader:
    '''
    A class to handel the JSON output from
        http://arohatgi.info/WebPlotDigitizer/
    It is taken from here
        https://github.com/ankitrohatgi/wpd-python
    '''

    def __init__(self, filename):
        with open(filename) as data_file:
            self.data = json.load(data_file)

    def getDatasetCount(self):
        return len(self.data['wpd']['dataSeries'])

    def getDatasetByIndex(self, index):
        return self.data['wpd']['dataSeries'][index]

    def getDatasetByName(self, name):
        return [x for x in self.data['wpd']['dataSeries']
                if x['name'] == name][0]

    def getDatasetNames(self):
        return [x['name'] for x in self.data['wpd']['dataSeries']]

    def getDatasetValues(self, dataset):
        values = []
        for val in dataset['data']:
            values.append(val['value'])
        return np.array(values)
