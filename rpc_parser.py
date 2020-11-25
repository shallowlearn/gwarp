'''
This work was supported by the Intelligence Advanced
Research Projects Activity (IARPA) via Department of
Interior / Interior Business Center (DOI/IBC) contract
number D17PC00280. The U.S. Government is authorized to
reproduce and distribute reprints for Governmental purposes
notwithstanding any copyright annotation
thereon. Disclaimer: The views and conclusions contained
herein are those of the authors and should not be
interpreted as necessarily representing the official
policies or endorsements, either expressed or implied, of
IARPA, DOI/IBC, or the U.S. Government.

Author: Bharath Comandur, cjrbharath@gmail.com
Date: 11/24/2020

Modified code written by Dr. Tanmay Prakash

'''

import re
import numpy as np


RPB_KEYS = [
    ('errBias','errBias'),
    ('errRand','errRand'),
    ('LINE_OFF','lineOffset'),
    ('SAMP_OFF','sampOffset'),
    ('LAT_OFF','latOffset'),
    ('LONG_OFF','longOffset'),
    ('HEIGHT_OFF','heightOffset'),
    ('LINE_SCALE','lineScale'),
    ('SAMP_SCALE','sampScale'),
    ('LAT_SCALE','latScale'),
    ('LONG_SCALE','longScale'),
    ('HEIGHT_SCALE','heightScale'),
    ('LINE_NUM_COEFF','lineNumCoef'),
    ('LINE_DEN_COEFF','lineDenCoef'),
    ('SAMP_NUM_COEFF','sampNumCoef'),
    ('SAMP_DEN_COEFF','sampDenCoef')]

RPB_EXTRA_KEYS = ['satId', 'bandId', 'SpecId', 'errBias', 'errRand']


# Load the parameters from an RPB file into a dict
def read_file_into_dict(rpcfile):
    f = open(rpcfile)
    d = {}
    flagGroup = False
    line = ''
    groupName = None
    for line_ in f:
        # Remove all whitespace
        line += line_.strip().replace(' ','').replace('\t','')
        if (not line.startswith('BEGIN_GROUP')
            and (not line.startswith('END_GROUP'))
            and (line.find(';') < 0)):
            # Keep appending until we have a ;
            continue
        if line.startswith('END;'):
            break
        # In case the line has stuff beyond the ;
        lines = line.split(';',1)
        if len(lines) > 1:
            line, lineNext = lines
        else:
            line = lines[0]
            lineNext = ''
        var, val = line.replace(';','').split('=')
        if flagGroup:
            if var == 'END_GROUP':
                assert groupName == val
                flagGroup = False
            else:
                d[groupName][var] = val
        else:
            if var == 'BEGIN_GROUP':
                groupName = val
                d[groupName] = {}
                flagGroup = True
            else:
                d[var] = val
        line = lineNext
    f.close()
    return d


def load_rpb(rpcfile):
    
    d = read_file_into_dict(rpcfile)
    rpc_coefs = {}
    p = re.compile('\(' + '([^,]+),'*19 + '([^,]+),*\)')
    for key, rpb_key in RPB_KEYS:
        if key.endswith('COEFF'):
            poly_coefs = [float(x) 
                          for x in p.match(d['IMAGE'][rpb_key]).groups()]
            poly_coefs = np.array(poly_coefs)
            rpc_coefs[key] = poly_coefs
        else:
            rpc_coefs[key] = float(d['IMAGE'][rpb_key])
    rpc_coefs['satId'] = d['satId']
    rpc_coefs['SpecId'] = d['SpecId']
    rpc_coefs['bandId'] = d['bandId']

    return rpc_coefs
