####################
# Written By Mattias Juhl
# Things to improve:
#  Need toadd a save name file on the export button.
# remove button inline, this really shouldn't be too hard. To be placed underneath the load button
# Can change plots
#######################


# Bugs fixed in this version
############
#       Saving bug, only saving one file
#       Reflection not working
#       Saving adding a tab at the end of the save file

# Changes from 2.4:
# 3
#       Added new input and output, based on a dictionary. This is in a file called ImportQSSFiles.
#       Labels, if can'tfind them will list number
##########


# Changes from 2.3:
# 3
#       B is now correct, updated using email from TT (2014)
#       Output now includes scalled data
#       Generation is no calculated differently for both PC and PL, previously it used Delta n from PC in both cases.
#       Graphs now auto scale after loggling
#       Legend button added for easier workings
##########


from numpy import *
from models.ConstantsClass import *
from CanvasClass import *


class QSSPL_Analyser(wx.App):

    """docstring for PVapp"""

    def OnInit(self):
        self.legacy_view = GUIController(None)
        self.legacy_view.Show()

        self.controller = PlaceholderController(self.legacy_view)

        self.data_panel_controller = DataPanelHandler(
            self.legacy_view.Data,
            self.legacy_view.data_panel,
            self.legacy_view
        )
        return True
