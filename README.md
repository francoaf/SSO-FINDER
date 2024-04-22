# SSO-FINDER
Pipeline for the detection of Solar System Objects in astronomical images.
SSO-FINDER finds a catalogue of moving objects (SSOs) by searching for sources in a list of catalogues, having the corresponding reference times (in MJD), that move from catalogue to catalogue.
A source is flagged as potential candidate if its proper motion and position angle is the same within given errors.

The file sso-finder.py contains the SSO class in which the "getit" function is the one useful for the SSO recovery.  
The function call is:

```
outralist, outdelist, outcatlist, times, CandidateSSO = getit(logger, dir, FileOFCat, FileOFRefTime, minDist, maxDist, ErrPosAngle, ErrPropMoti, format="txt", racolumn=0, decolumn=1, MinimumNumberCATON=1, MinimumCat=3, info=True)
```

where in output you have:
* _outralist_:     Right Ascention of candidates
*_outdelist_:     Declination of candidates
* _outcatlist_:    List of catalogues in which candidates are found
* _times_:         List of times corresponding to outcatlist
* _CandidateSSO_:  Candidates list (indexes in outcatlist)

An articole describing the usage and the capabilities of the SSO-FINDER is now in preparation.

Licensed under version 3 of the GNU Affero General Public License.
