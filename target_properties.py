#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import astropy.units as u
from astropy.table import Table, QTable, Column
from astropy import table
from astropy.coordinates import SkyCoord, Distance, Latitude, Longitude, Distance
from astropy.time import Time

from astroquery.simbad import Simbad


# In[ ]:


names = Table.read('targetnames.txt', format='ascii.no_header', delimiter='|', names=['target'])


# In[ ]:


TargetSimbad = Simbad()
TargetSimbad.add_votable_fields('propermotions', 'distance', 'parallax', 'flux(B)', 'flux(V)', 'rv_value', 'sp')
TargetSimbad.add_votable_fields('ra(2;A;ICRS;J2000;2000)', 'dec(2;D;ICRS;J2000;2000)')


# In[ ]:


targets_download = TargetSimbad.query_objects(names['target'])


# In[ ]:


targetsin = table.hstack([names, targets_download], join_type='exact')
# translate units to astropy
for col in targetsin.colnames:
    if str(targetsin[col].unit) == '"h:m:s"':
        targetsin[col] = Longitude(targetsin[col], unit=u.hourangle)
    elif str(targetsin[col].unit) == '"d:m:s"':
        targetsin[col] = Latitude(targetsin[col], unit=u.deg)
# For Distance_distance, units are in separate column 


# In[ ]:


targets = QTable(targetsin)


# In[ ]:


targets['Distance_unit'][targets['Distance_unit'] == ''] = 'pc'
targets['distance'] = u.Quantity([row['Distance_distance'] * u.Unit(row['Distance_unit']) for row in targets])
# Some distances suddenly become 0, where the input was masked. We deal with those next.


# Distances are typically taken from GAIA, but some bright stars are missing. 
# For those, we calculate the distance from the parallax (from HIPPARCOS)
ind = targets['Distance_distance'].mask
ind = targets['Distance_perr'].mask
targets['distance'][ind] = targets['PLX_VALUE'].to(u.pc, equivalencies=u.parallax())[ind]
targets['Distance_method'][ind] = 'paral'
targets['Distance_perr'][ind] = ((targets['PLX_VALUE'] - targets['PLX_ERROR']).to(u.pc, equivalencies=u.parallax()) 
                                 - targets['distance'])[ind]
targets['Distance_merr'][ind] = ((targets['PLX_VALUE'] + targets['PLX_ERROR']).to(u.pc, equivalencies=u.parallax()) 
                                 - targets['distance'])[ind]
targets['Distance_bibcode'][ind] = targets['PLX_BIBCODE'][ind]
# finally: Clean up non-quantity input columns
targets.remove_columns(['Distance_distance', 'Distance_unit'])


# In[ ]:


# Convert byte objects into strings, because we can't save byte objetcs to fits later
for c in targets.colnames:
    if targets[c].dtype == object:
        if type(targets[c][0]) is str:
            targets[c] = [s for s in targets[c]]
        else:
            targets[c] = [s.decode() for s in targets[c]]


# In[ ]:


# Because of https://github.com/astropy/astropy/issues/13041
# input must be a Quantity, not a maskedQuantity.
# So, we fill in 0 otherwise (no proper motion).
# Currently, all numbers are valid though, so the "filled" does not change anything,
# it just converts form MaskedQuantity to Quantity.
targets['coord'] = SkyCoord(ra=targets['RA_2_A_ICRS_J2000_2000'],
                            dec=targets['DEC_2_D_ICRS_J2000_2000'],
                            distance=targets['distance'],
                            pm_ra_cosdec=targets['PMRA'].filled(0),
                            pm_dec=targets['PMDEC'].filled(0),
                            frame='icrs', obstime='J2000', equinox='J2000')


# In[ ]:


targets.write('targets.csv', format='ascii.ecsv', overwrite=True)

