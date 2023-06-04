'''
Utilities
Add things here as needed
'''

# Standard imports
import sciris as sc
import pandas as pd
import numpy as np


def set_font(size=None, font='Libertinus Sans'):
    ''' Set a custom font '''
    sc.fonts(add=sc.thisdir(aspath=True) / 'assets' / 'LibertinusSans-Regular.otf')
    sc.options(font=font, fontsize=size)
    return
