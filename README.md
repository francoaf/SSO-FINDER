# SSO-FINDER
Pipeline for the detection of Solar System Objects in astronomical images.
SSO-FINDER finds a catalogue of moving objects (SSOs) by searching for sources in a list of catalogues, having the corresponding reference times (in MJD), that move from catalogue to catalogue.  
An object is flagged as potential candidate if its proper motion and position angle is the same within given errors.

The file sso-finder.py contains the SSO class in which the "sso_finder" function is the one useful for the SSO recovery.  
The function call is:

```
outralist, outdelist, outcatlist, times, CandidateSSO = sso_finder(logger, dir, FileOFCat, FileOFRefTime, minDist, maxDist, ErrPosAngle, ErrPropMoti, format="txt", racolumn=0, decolumn=1, MinimumNumberCATON=1, MinimumCat=3, info=True)
```

where in output you have:
* _outralist_:     Right Ascention of candidates
*_outdelist_:     Declination of candidates
* _outcatlist_:    List of catalogues in which candidates are found
* _times_:         List of times corresponding to outcatlist
* _CandidateSSO_:  Candidates list (indexes in outcatlist)

An articole describing the usage and the capabilities of the SSO-FINDER is now in preparation.

### Licensed under version 3 of the GNU Affero General Public License.

This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 3.0 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
