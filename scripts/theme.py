"""Reusable publication theme; physical sizes are in millimetres."""
import matplotlib as mpl
from matplotlib import font_manager
AG = '#527D91'
BORZOI = '#747D82'
REFERENCE = '#414B50'
INSERTION = '#5D8884'
LIGHT = '#CBD1D3'
TEXT = '#252C30'
COLORS = {'AlphaGenome':AG, 'Borzoi':BORZOI, 'Reference':REFERENCE, 'Off-reference':INSERTION}
MARKERS = {'AlphaGenome':'o', 'Borzoi':'s', 'Reference':'o', 'Off-reference':'D'}
WIDTH_MM = 170
def apply():
    try:
        font_manager.findfont('Arial', fallback_to_default=False)
        family = 'Arial'
    except ValueError:
        family = 'DejaVu Sans'
    mpl.rcParams.update({'font.family':family,'font.size':7.5,'axes.labelsize':8,
        'axes.titlesize':8.5,'xtick.labelsize':7.2,'ytick.labelsize':7.2,
        'legend.fontsize':7.2,'axes.linewidth':.8,'lines.linewidth':.9,
        'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False,
        'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white',
        'text.color':TEXT,'axes.labelcolor':TEXT,'xtick.color':TEXT,'ytick.color':TEXT,
        'svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42,'text.usetex':False})
