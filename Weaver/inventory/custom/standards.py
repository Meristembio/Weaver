from Bio.Restriction.Restriction_Dictionary import rest_dict
from ..models import RestrictionEnzyme

sapi_enzyme_name = 'sapi'
bsmbi_enzyme_name = 'bsmbi'
sapi = None
bsmbi = None

if sapi_enzyme_name in rest_dict:
    sapi = RestrictionEnzyme(name=sapi_enzyme_name)
if bsmbi_enzyme_name in rest_dict:
    bsmbi = RestrictionEnzyme(name=bsmbi_enzyme_name)

ligation_standards = {
    'loop': {
        'name': 'Loop',
        'family': 'golden_gate',
        'domestication': {
            'enzyme': sapi,
            'receiver': {
                'name': 'pL0R-mRFP1',
                'ohs': {
                    'oh5': 'TCC',
                    'oh3': 'CGA'
                },
            },
            'enzymes': [
                'aari',
                'bsai',
                'sapi'
            ],
        },
        'ohs': {
            'l0': {
                'a': {
                    'name': 'A',
                    'oh': 'GGAG'
                },
                'b': {
                    'name': 'B',
                    'oh': 'TACT'
                },
                'c': {
                    'name': 'C',
                    'oh': 'AATG'
                },
                'd': {
                    'name': 'D',
                    'oh': 'AGGT',
                    'tc': True
                },
                'e': {
                    'name': 'E',
                    'oh': 'GCTT',
                    'stop': True
                },
                'f': {
                    'name': 'F',
                    'oh': 'CGCT'
                },
            },
            'l1': {
                'alpha': {
                    'name': 'GAMMA',
                    'oh': 'ATG'
                },
                'beta': {
                    'name': 'BETA',
                    'oh': 'GCA'
                },
                'gamma': {
                    'name': 'GAMMA',
                    'oh': 'TAC'
                },
                'epsilon': {
                    'name': 'EPSILON',
                    'oh': 'CAG'
                },
                'omega': {
                    'name': 'OMEGA',
                    'oh': 'GGT'
                },
            },
        },
    },
    'gb': {
        'name': 'Golden Braid 2.0',
        'family': 'golden_gate',
        'domestication': {
            'enzyme': bsmbi,
            'receiver': {
                'name': 'pUPD',
                'ohs': {
                    'oh5': 'CTCG',
                    'oh3': 'TGAG'
                },
            },
            'enzymes': [
                'bsai',
                'btgzi',
                'bsmbi'
            ],
        },
        'ohs': {
            'l0': {
              'a1_5': {
                'name': 'A1 5\'',
                'oh': 'GGAG'
              },
              'a2_5': {
                'name': 'A2 5\'',
                'oh': 'TGAC'
              },
              'a3_5': {
                'name': 'A3 5\'',
                'oh': 'TCCC'
              },
              'b1_5': {
                'name': 'B1 5\'',
                'oh': 'TACT'
              },
              'b2_5': {
                'name': 'B2 5\'',
                'oh': 'CCAT'
              },
              'b3_5': {
                'name': 'B3 5\'',
                'oh': 'AATG'
              },
              'b4_5': {
                'name': 'B4 5\'',
                'oh': 'AGCC',
                'tc': True
              },
              'b5_5': {
                'name': 'B5 5\'',
                'oh': 'TTCG',
                'tc': True
              },
              'b6_5': {
                'name': 'B6 5\'',
                'oh': 'GCTT',
                'stop': True
              },
              'c1_5': {
                'name': 'C1 5\'',
                'oh': 'GGTA'
              },
              'c1_3': {
                'name': 'C1 3\'',
                'oh': 'CGCT'
              },
            },
            'l1': {
            },
        },
    },
}
