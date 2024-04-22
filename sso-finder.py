#
# Copyright (C) 2012-2024 Nucita A.A., Franco A., Sacquegna S.
# Department of Mathematics and Physics ``E. De Giorgi'' , University of Salento, Via per Arnesano, CP-I93, I-73100, Lecce, Italy
# INFN, Sezione di Lecce, Via per Arnesano, CP-193, I-73100, Lecce, Italy
# INAF, Sezione di Lecce, Via per Arnesano, CP-193, I-73100, Lecce, Italy
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3.0 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#
# File: SSOFinder.py
# Fate: 15/06/22
# Authors: Nucita A.A., Franco A., Sacquegna S.
# Note:  A lass to find MOVING objects in ALREADY calibrated catalogue
#
# INPUT:
#     # A logger objects
#     logger: logger object
#     # The input folder
#     dir: string, path to the data
#     # PIVOT catalogue
#     PivotCatalogue: string, filename of existing catalogue
#     # Reference time associated to PIVOT catalogue
#     PivotReferenceTime: float, MJD, days
#     # Input catalogue list
#     InputListCatalogue: list of string, filenames of existing catalogues
#     # Reference time list associated to input catalogue list (same number of elements as InputListCatalogue)
#     InputListTimes: list of float, MJD, days (same number of elements as InputListCatalogue).
#     # Minimum Distance (arcsec) for skipping SSO candidateship
#     minDist: float, minimum distance for matching  (arcsec)
#     # Maximum distance (arcsec) to searh for each PIVOT target
#     maxDist: maxim distance for matching (arcsec)
#     # Error on position angle (degrees)
#     ErrPosAngle: float, error in angle, degrees
#     # Error on proper motion (arcsec/minutes)
#     ErrPropMoti: float, error in proper motion, degrees
#     # Minimum number of catalogues that MUST contain an object in order to declare it NOT a good candidate
#     MinimumNumberCATON: integer.
#     # format: string. Must be "txt" or "fits" depending on the inpyt data. If ASCII files are provided, data must be
#               formatted with two columns having the meaning of right ascension (in degrees) and declination
#               (in degrees). Id fits files are passed, these file must contain at least two columns RA and DEC.
#     # extension number in case a
#     ext
#
# OUTPUT: getit functions returns:
#         rap, dep: coordinates numpy arrays of objects in PIVOT catalogue
#         ralist, delist: list of coordinates of numpy arrays. Each element of the list is a set of coordinates
#                         corresponding to one of the input catalogues.
#         return -1, -1, -1,  if a severe error is raised
# -------------------------------------------------------------------------------

import numpy
import os
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
import scipy
import matplotlib.pyplot as plt
from pdb import set_trace as stop
import pandas as pd

pd.options.mode.chained_assignment = None


class SSO(object):
    """

    """

    def __init__(self):
        """


        """
    # ----------------------------------------------------------------------------
    def getit(self, logger, dir, FileOFCat, FileOFRefTime, minDist, maxDist, ErrPosAngle, ErrPropMoti, format="txt",
              racolumn=0, decolumn=1, MinimumNumberCATON=1, MinimumCat=3, info=True):
        '''

        FIND a catalogue of moving objects (SSO) by searching for sources in FileOFCat (list of filename of catalogues) having times FileOFRefTime (in MJD) that move from catalogue to catalogue.
        A source is associated if the speed and position angle is the same within given errors.


        For explanation on the input and usage, see the comment section in the header of this file.
        :param logger
        :param dir directory where the catalogues are stored
        :param FileOFCat: list fo filneame of catalogues. First entry is assumed to be the pivot
        :param FileOFRefTime: list of reference times. On time (in days) for each cataloge in FileOfCat. First entry is assumed to be the reference time of the pivot
        :param minDist: minimum distance (i.e. the associated astrometric error between catalogues) to consider sources in different catalogues as the same sources.
        :param maxDist: maximum distance (for each pivot source) that the code consider for the association of sources. It is the maximum distance travelled by the pivot (if a real moving object) in the time between the first and last catalogue
        :param ErrPosAngle: error in dgerees according to which two sources are considered aligned
        :param ErrPropMoti: error in proper motion (percentage) according to which to sources are considered associated
        :param MinimumNumberCATON: 1 by default. In the present version of the code the input catalogues are firs purged in order to remove all teh sources that are within minDist.
        :param format txt by default. No check is done on varius input. There is the plan to extend the reading routine also to sextractor ldac fits files
        :param info to print a few information on screen during the run
        :parm  racolum integer corresponding to the ra column to be read from the ASCII file (in degree)
        :parm  decolum integer corresponding to the dec column to be read from the ASCII file (in degree)
        :param MinimumCat minimum number of catalogues in whcih a source was identified as a moving object in order to be considered a good candidate. Default 3
        :return:
        '''
        # selfs
        self.logger = logger
        self.dir = dir
        self.FileOFCat = FileOFCat
        self.FileOFRefTime = FileOFRefTime
        self.minDist = minDist
        self.maxDist = maxDist
        self.ErrPosAngle = ErrPosAngle
        self.ErrPropMoti = ErrPropMoti
        self.format = format
        self.MinimumNumberCATON = MinimumNumberCATON
        self.MinimumCat = MinimumCat
        self.info = info
        self.racolumn = racolumn
        self.decolumn = decolumn

        logger.info('getit:> reading input list files...')
        p = os.path.join(self.dir, self.FileOFCat)
        if (not os.path.isfile(p)):
            logger.info('getit:> File does not exist. Severe stop. Check: %s', str(p))
            return -1, -1, -1, -1
        else:
            tmpcats = numpy.loadtxt(p, dtype=str, usecols=(0), unpack=True, comments='#')

        p = os.path.join(self.dir, self.FileOFRefTime)
        if (not os.path.isfile(p)):
            logger.info('getit:> File does not exist. Severe stop. Check: %s', str(p))
            return -1, -1, -1, -1
        else:
            tmpref = numpy.loadtxt(p, dtype=float, usecols=(0), unpack=True, comments='#')

            # All catalogues filename
            InputListALLCatalogue = numpy.ndarray.tolist(tmpcats)

            # Number of catalogues (including pivot)
            NAllCatalogue = len(InputListALLCatalogue)

            # Pivot catalogue filename
            PivotCatalogue = tmpcats[0]

            # Name of N-1 catalogues (pivot excluded)
            InputListCatalogue = numpy.ndarray.tolist(tmpcats[1:])

            # Number of N-1 catalogues (pivot excluded)
            NCatalogues = len(InputListCatalogue)

            # Pivot Reference time in days
            PivotReferenceTime = tmpref[0]

            # Other catalogues reference time in days
            InputListTimes = tmpref[1:]


            Ntimes = len(InputListTimes)

        # Prepare an array containing the distance in time (in hour) between each catalogue and teh pivot catalogue
        ImageTimes = tmpref

        # DeltaT MJD
        DeltaT = (InputListTimes - PivotReferenceTime) * 24.0e0

        # Maximum separation to be searched for
        max_sep = self.maxDist * u.arcsec
        # Minimum separation to be searched for. two stars below this separation im
        min_sep = self.minDist * u.arcsec

        logger.info('getit:> reading ASCII files...')

        # List of arrays (coordinates for all catalogues)
        ralistAll = []
        delistAll = []

        # At the end of the code, it is necessary to recover all the useful information for each sso found (as magnitude and other parameters found by Sextractor).
        # Instead of moving forward all these information, we simply move the name of the catalogue from which that SSO (in that particular instant of time) was found.
        # Therefore, we build a list of lists of Catalogue filenames. catlistALL contains N lists of names of catalogues. N is the number of input catalogues.
        # The first element of the list associated to catalogue 1) is itself a list containing M times the name of catalogue 1. M is the number of objects found by sextractor.

        catlistALL = []

        for ilist in InputListALLCatalogue:
            c = os.path.join(self.dir, ilist)
            if (not os.path.isfile(c)):
                logger.info('getit:> File does not exist. Severe stop. Check: %s', str(c))
                return -1, -1, -1, -1
            # rac, dec = numpy.loadtxt(c, usecols=(0, 1), unpack=True)
            rac, dec = numpy.loadtxt(c, usecols=(self.racolumn, self.decolumn), unpack=True)
            self.makeds9reg(logger, os.path.join(self.dir, c.split('.txt')[0] + '.reg'), rac, dec, sys='fk5',
                            colour='red')

            ralistAll.append(rac)
            delistAll.append(dec)

            # For this catalogue, creating a list containing M times the name of the catalogue itself (ilist) and appending this list to catlistALL
            repeatedCatNameList = [ilist] * len(rac)
            # Append this list of repeated names to the catlistALL
            catlistALL.append(repeatedCatNameList)

        # Final SSO candidate list. It has lists of 4 elements each with the idex of the matched source moviving in each catalogue. Firts entry is for the pivot. Other entries for the subsequent
        # catalogues. -1 means that the sources was not found in that catalogue
        CandidateSSO = []

        # Deleting sources within MinDist
        logger.info('getit:> removing sources within min_sep...')

        index_to_be_removed = self.cross_correlate_multi_cat(NAllCatalogue, ralistAll, delistAll, sep=min_sep)

        # Purging input catalogues
        for icat in range(0, len(ralistAll)):
            ralistAll[icat] = numpy.delete(ralistAll[icat], index_to_be_removed[icat])
            delistAll[icat] = numpy.delete(delistAll[icat], index_to_be_removed[icat])
            catlistALL[icat] = numpy.delete(catlistALL[icat], index_to_be_removed[icat])

        # Considering pivot
        rap = ralistAll[0]
        dep = delistAll[0]

        cap = catlistALL[0]

        # this stores the pivot index in the pivotindex
        pivotindex = numpy.arange(0, len(rap))

        # All catalogues, except the pivot. Each element of the list is a catalogue
        ralist = ralistAll[1:]
        delist = delistAll[1:]

        catlist = catlistALL[1:]

        # write ds9 region for each catalogue
        # Pivot catalogue

        # Check again that no star in pivot has a corresponding associated entry (within min_sep) in the other catalogues)
        # Second method: using separation and then masking the output with res = (sep.arcsecond) * u.arcsec < min_sep. Then counting the len of reslist. Repeating for each catalogues, increase oncat if an association is found

        logger.info('getit:> searching candidates with the input constraints...')
        for indexpivot, ra, dec in zip(pivotindex, rap, dep):

            cpivot = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')

            # Run over the catalogues and search for cross correlated sources within min_sep
            # This check is now useless since all the stars checked across the catalogues do not correlate within min_sep. We have already purged all the catalogues in this respect
            # In the old code, the check was performed at this level and the variable oncat increased if a source is found on more that MinimumNumberCATON catalogues. The source is then excluded
            # Here, sìince we already did the check, oncat is always set to 0
            oncat = 0

            if (oncat >= self.MinimumNumberCATON):
                # The pivot index is not considered as a possible target
                if (self.info):
                    print("PIVOT is not SSO: ", str(indexpivot))
                continue
            else:
                # The pivot index is considered as a possible target and is appended to a tmp list containing its index in the pivot catalogue
                if (self.info):
                    print('')
                    print("PIVOT is flagged to be a potential SSO: ", str(indexpivot))

                # prepare the tmp list and append indexpivot to it if oncat = 0
                tmp = []
                tmp_new = []
                if (self.info):
                    print('PIVOT index: ', indexpivot)
                    print('RA (deg), DEC (deg): ', rap[indexpivot], dep[indexpivot])
                tmp.append(indexpivot)
                tmp_new.append(indexpivot)

                # For this pivot and indexpivot, Let's go across the catalogues and search for pairs
                listPairProperMotion = []
                ListPair = []
                ListPairIndices = []
                ListPairPositionAngle = []

                for i in range(NCatalogues):
                    ci = SkyCoord(ra=ralist[i] * u.deg, dec=delist[i] * u.deg, frame='icrs')
                    sep = cpivot.separation(ci)

                    # get sources only within distMax (maximum velocity measurable)
                    # get sources only candidates with distance larger than minDist, thus avoiding to consider pivot (if any) in other catalogues.

                    res = ((sep.arcsecond) * u.arcsec <= max_sep) & ((sep.arcsecond) * u.arcsec >= min_sep)
                    # Get the coordinates of the matched sources
                    ci_matches = ci[res]
                    if (self.info):
                        print('')
                        print('pivot: ', indexpivot, 'catalogue: ', i, '# associated to pivot: ', len(ci_matches))
                        for imat, ci_mat in enumerate(ci_matches):
                            print('match :', imat, 'RA: ', ci_mat.ra.deg, 'DEC: ', ci_mat.dec.deg)

                    # Append to list of pairs
                    ListPair.append(ci_matches)

                    # convert bolean list of associations to array of indices (contains indices of the catalogue that are associated to the pivot)
                    ListPairIndices.append([i for i, x in enumerate(res) if x])

                    # Evaluate distance in arcseconds of pair
                    di_matches = sep.arcsecond[res]

                    # and evaluate the proper motion arcsec/hour and append to list of proper motion
                    mi_matches = di_matches / DeltaT[i]

                    # append the mi_matches to the listPairProperMotion
                    listPairProperMotion.append(mi_matches)

                    # Evaluate pair position angle
                    # evaluate delta Dec = dec catalogue object - dec pivot
                    deltaDec = ci_matches.dec.degree - cpivot.dec.degree
                    # evaluate delta Ra = ra catalogue object - ra pivot
                    deltaRa = ci_matches.ra.degree - cpivot.ra.degree
                    # thetai = arctan(deltaDec/deltaRa)
                    thetai_matches = numpy.arctan2(deltaDec, deltaRa)

                    # get Omegai (if thetai < 0  omegai = thetai+2*pi), if thetai >0, omegai=thetai
                    omegai_matches = numpy.zeros(len(thetai_matches))
                    k = 0
                    for th in thetai_matches:
                        if th < 0:
                            om = th + 2 * numpy.pi
                        else:
                            om = th
                        omegai_matches[k] = om
                        k = k + 1

                    # Append to ListPairPositionAngle
                    # ListPairPositionAngle.append(thetai_matches)
                    ListPairPositionAngle.append(omegai_matches)

                # At this stage, for each pivot object, we have a list of sources around it with distance (wrt the pivot)
                # in the annulus with radii distMIN, and distMAX.
                # we now find the objects in the separate catalogues that have the same mui and omegai
                # get indices in catalogues

                if (self.info):
                    print('getit:> associating...')
                array_tmp = numpy.zeros((NCatalogues, NCatalogues, NCatalogues), dtype=int)
                array_tmp.fill(-1)

                for ifirst_cat in range(0, NCatalogues - 1):
                    for jsecond_cat in range(ifirst_cat + 1, NCatalogues):
                        # setting variables icat1 and icat2 just for a dummy match with the old loop
                        icat1 = ifirst_cat
                        icat2 = jsecond_cat
                        # setting the maximum values of the speed and angle to match. These values are used to search for sources with speed and position angle within these limits
                        maxmu = self.ErrPropMoti
                        maxang = self.ErrPosAngle
                        # setting the temporary maximum values of the speed and angle to match. These values are used to search for the source (among those with speed and position angle within the above limit) with minimum requirements
                        maxmu_tmp = maxmu
                        maxang_tmp = maxang
                        # Starting all the combinations
                        if (self.info):
                            print('correlating speed and position angle catalogues: ', icat1, icat2)
                        # By default (as a firts approximation) a candidate SSO is flagged as NOT BEING an SSO, i.e. indcat1 and indact2 are set to -1. If two sources in cat 1 and cat2
                        # are within the speed and angle limits, the variables are updated with the i and j indices asscoiated to those sources in teh respective catalogues. i is for catalogue 1, j is for catalogue 2
                        indcat1 = -1
                        indcat2 = -1
                        for i in range(len(listPairProperMotion[icat1])):
                            for j in range(len(listPairProperMotion[icat2])):
                                # percentage error in velocity
                                # dist = numpy.abs(listPairProperMotion[icat1][i] - listPairProperMotion[icat2][j]) / (numpy.abs(max(listPairProperMotion[icat1][i], listPairProperMotion[icat2][j])))
                                dist = numpy.abs(
                                    listPairProperMotion[icat1][i] - listPairProperMotion[icat2][j]) / numpy.abs(
                                    listPairProperMotion[icat1][i])
                                distang = numpy.abs(
                                    ListPairPositionAngle[icat1][i] - ListPairPositionAngle[icat2][j]) * 180. / numpy.pi

                                if (dist <= maxmu and distang <= maxang):
                                    if (dist <= maxmu_tmp and distang <= maxang_tmp):
                                        maxmu_tmp = dist
                                        maxang_tmp = distang
                                        indcat1 = i
                                        indcat2 = j
                        if (self.info):
                            print("associated indices:", indcat1, indcat2)

                        #############################
                        # If indexes indcat1 and indcat2 turn to be different from -1 then there is an association
                        if (indcat1 >= 0 and indcat2 >= 0):
                            array_tmp[icat1, icat2, icat1] = ListPairIndices[icat1][indcat1]
                            array_tmp[icat1, icat2, icat2] = ListPairIndices[icat2][indcat2]
                            if (self.info):
                                print(ListPairIndices[icat1][indcat1], listPairProperMotion[icat1][indcat1],
                                      ListPair[icat1][indcat1].ra.degree, ListPair[icat1][indcat1].dec.degree)
                            tmp.append(ListPairIndices[icat1][indcat1])
                            if (self.info):
                                print(ListPairIndices[icat2][indcat2], listPairProperMotion[icat2][indcat2],
                                      ListPair[icat2][indcat2].ra.degree, ListPair[icat2][indcat2].dec.degree)
                            tmp.append(ListPairIndices[icat2][indcat2])
                        else:
                            array_tmp[icat1, icat2, icat1] = -1
                            array_tmp[icat1, icat2, icat2] = -1
                            tmp.append(-1)
                            tmp.append(-1)

                # At this level the association is done and stored in the array array_tmp. Now we recover the array_tm
                # association and store into a list
                # run over the catalogues indices 0, -Ncatalogue -1
                # the index of the associated star in the catalogue icatalogue_back is stored in the array
                # array_tmp[:,:,icatalogue_back].
                # This is a Ncatalogue*Ncatalogue matrix with all elements = -1 if no association is found for that
                # catalogue, and (hopefully) with a few components different from -1 if an association is found
                # In this case, the integer stored (and different from -1) is the index of the associated star in the
                # catalogue # icatalogue_back!
                # We perform to checks: for each icatalogue_back all the elements = -1 ---> no association and -1 is
                # appended to the tmp_new list
                # If a few components are different form -1, We check that these componente store the same index.
                # IT MUST BE THE SAME. IF SO the integer is stored in the tmp_new list, otherwise an error is raised,
                # meaning that the association failed. In this case we append -1

                for icatalogue_back in range(0, NCatalogues):
                    # get the matrix for this catalogue
                    icatalogue_matrix = array_tmp[:, :, icatalogue_back]
                    # check for -1
                    res_icatalogue_matrix = numpy.where(icatalogue_matrix != -1)
                    # extract the matrix. It is a  zero element list if all the elements of the matrix associate to
                    # the catalogue are -1. Otherwise, the len is differert from 0. In this case, the result isn a
                    # numpy array containing the ndices of the correlates sourec in the particular catalogue
                    extract_icatalogue_matrix = icatalogue_matrix[res_icatalogue_matrix]
                    # if len = 0 all the elements are -1 and we append to tmp_new -1
                    if len(extract_icatalogue_matrix) == 0:
                        tmp_new.append(-1)
                    else:
                        # check if all elemnts are equal. If equal, append the first elements
                        if (self.all_equal(extract_icatalogue_matrix.tolist())):
                            tmp_new.append(extract_icatalogue_matrix[0])
                        else:
                            # print("ERROR: different indices in associating catalogues: ", icatalogue_back,'for pivot: ', indexpivot, extract_icatalogue_matrix)
                            tmp_new.append(-1)

                # Appending candidate (called "tmp") to the CandidateSSO list
                tmp = tmp_new
                if (self.info):
                    print('list associated indices: ', tmp)
                CandidateSSO.append(tmp)

        ################################################################
        # Keep only candidates with >= MinimumCat elements not negative (-1)
        ###############################################################
        logger.info('getit:> Purging candidate SSO with MinimumCat...')
        CandidateSSOtmp = []
        for tmp in CandidateSSO:
            count = 0
            for i in range(NCatalogues + 1):
                if tmp[i] >= 0:
                    count = count + 1
            if (count >= (self.MinimumCat)):
                CandidateSSOtmp.append(tmp)
        CandidateSSO = CandidateSSOtmp

        ################################################################
        # Rearrange output lists
        ###############################################################
        logger.info('getit:> Appending...')
        outcatlist = []
        outralist = []
        outdelist = []
        outcatlist = []
        outcatlist.append(cap)
        for i in catlist:
            outcatlist.append(i)
        outralist.append(rap)
        for i in ralist:
            outralist.append(i)
        outdelist.append(dep)
        for i in delist:
            outdelist.append(i)
        # Times.. rearrange times as a numpy array
        times = []
        times.append(PivotReferenceTime)
        for i in InputListTimes:
            times.append(i)
        times = numpy.asarray(times)

        logger.info('getit:> --------------------------------------')
        logger.info('%s SSO candidates have been found', str(len(CandidateSSO)))
        logger.info('getit:> --------------------------------------')

        return outralist, outdelist, outcatlist, times, CandidateSSO
