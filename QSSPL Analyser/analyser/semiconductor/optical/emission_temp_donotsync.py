
# this is an old class, that i'm not sure is read yet.

class Alpha_from_PL():

    """
    This class is for given a PL spectrum trying to
    determine the absorption coefficents
    """
    _cal_dts = dict(
        material='Si',
        temp=300.,  # temp in kelvin
        width=0.018,  # width in cm
        ni_author=None,  # author of intrinsic carrier density
        optics_k_author='Green2008',
        optics_n_author='Green2008',
        nxc=np.ones(10),  # the number of excess carrier with depth
        doping=1e16,  # the doping in cm^-3
        wafer_opitcs='polished',
        detection_side='front',
        )

    known_alpha = 1
    known_wavelength = 1
    wavelength_measured = np.array([1])
    PL = np.array([1])

    W = 0.018
    # Dictionaries
    wafer_optics_dic = {'polished': 'escprob_polished_schick1992',
                        'textured': 'escprob_textured_rudiger2007'}
    PL_Dection_side_depth = {'rear': 'Escape_rear',
                             'front': 'Escape_front'}
    wafer_opitcs = 'polished'
    detection_side = 'front'

    # this class should start by doing a point wise fit to the
    def PL_normtoBB(self):
        """
        Just a little function to remove black body stuff
        """

        BB = BlackBody.PhotonFlux(self.wavelength_measured, self._cal_dts['temp'])

        self.PLnBB = self.PL / BB

    def Guess_alpha(self):
        """
        Used to align the realtive PL measurement to calibrated alpha at
        a known wavelength. Then does a first guess at to what alpha is.
        """

        self.PL_normtoBB()
        norm_PL = np.interp(
            self.known_wavelength, self.wavelength_measured, self.PLnBB)

        # Find the proportionality in PL for the known alpha
        self.wavelength = np.array(
            (self.known_wavelength, self.known_wavelength))
        self.optics_abs_cof_bb = np.array((self.known_alpha, self.known_alpha))
        self.update_escape()

        # For checking the escape probability
        # plt.figure('test')
        # plt.plot(self.x, self.escapeprob[:,0])
        # plt.show()
        self.Calibrationconstant = norm_PL / self.known_alpha / \
            np.trapz(self.np * self.escapeprob[:, 0], self.x)

        self.PL_alpha = self.PLnBB / norm_PL * self.known_alpha

    def update_escape(self):
        """
        calculates the escape probability given alpha
        """

        getattr(self, self.wafer_optics_dic[self.wafer_opitcs])()

        self.escapeprob = getattr(
            self, self.PL_Dection_side_depth[self.detection_side])

    def Iterate_alpha(self, n=10):
        """
        A function (not verified) to determine determine alpha from a PL
        spectrum itterativly

        Does the following itteration n times:

        1. It assumes alpha, calcs the escape probaility, calculates PL
        2. The calc PL is compared to the real PL and alpha is updated


        a note:
        Kramers Kronig could be used here to provide a relationship for how
        alpha should behave
        "http://www.doria.fi/bitstream/handle/10024/96800/Conventional%20and%20nonconventional%20Kramers-Kronig%20analysis%20in%20optical%20spectroscopy.pdf?sequence=3"
        """

        for i in range(n):
            # print i
            self.wavelength = self.wavelength_measured
            self.optics_abs_cof_bb = self.PL_alpha

            self.update_escape()

            self.PL_alpha = self.PLnBB / self.Calibrationconstant /\
                np.trapz(self.np * self.escapeprob.T,
                         self.x,
                         axis=1)
