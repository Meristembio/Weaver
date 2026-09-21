const L0DPartsStandards = {
  ytk: {
    name: 'Yeast Toolkit',
    enzyme: 'bsmbi',
    bases_upto_snip: 'a',
    l0_receiver: {
        name: 'L0R-GFP',
        ohs: {
            'oh5': 'CTCC',
            'oh3': 'CGAG'
        },
    },
    domestication_enzymes:[
            'bsmbi',
            'bsai',
    ],
    default: {
      5: '12',
      3: '23'
    },
    ohs: {
        '81': {
            name: '8-1',
            oh: 'CCCT'
        },
        '12': {
            name: '1-2',
            oh: 'AACG'
        },
        '23': {
            name: '2-3',
            oh: 'TATG'
        },
        '3a3b': {
            name: '3a-3b',
            oh: 'TTCT',
            prev_bases: 'gg'
        },
        '34': {
            name: '3-4',
            oh: 'ATCC',
            stop: true
        },
        '4a4b': {
            name: '4a-4b',
            oh: 'TGGC'
        },
        '45': {
            name: '4-5',
            oh: 'GCTG'
        },
        '56': {
            name: '5-6',
            oh: 'TACA'
        },
        '67': {
            name: '6-7',
            oh: 'GAGT'
        },
        '78': {
            name: '7-8',
            oh: 'CCGA'
        },
        '8a8b': {
            name: '8a-8b',
            oh: 'CAAT'
        },
    },
  },
  loop: {
    name: 'Loop',
    enzyme: 'sapi',
    bases_upto_snip: 'a',
    l0_receiver: {
      name: 'pL0R-mRFP1',
      ohs: {
        oh5: 'TCC',
        oh3: 'CGA'
      },      
    },
    domestication_enzymes:[
        'aari',
        'bsai',
        'sapi'
    ],
    default: {
      5: 'a',
      3: 'b'
    },
    ohs: {
      a: {
        name: 'A',
        oh: 'GGAG'
      },
      b: {
        name: 'B',
        oh: 'TACT'
      },
      c: {
        name: 'C',
        oh: 'AATG'
      },
      d: {
        name: 'D',
        oh: 'AGGT',
        prev_bases: 'tc'
      },
      e: {
        name: 'E',
        oh: 'GCTT',
        stop: true
      },
      f: {
        name: 'F',
        oh: 'CGCT'
      },
    }
  },
  gb: {
    name:"GoldenBraid 2.0",
    enzyme: 'bsmbi',
    bases_upto_snip: 'a',
    domestication_enzymes:[
      'bsai',
      'bsmbi'
    ],
    l0_receiver: {
      name: 'pUPD',
      ohs: {
        oh5: 'CTCG',
        oh3: 'CGAG'
      },      
    },
    default: {
      5: 'a1_5',
      3: 'a2_5'
    },
    ohs: {
      a1_5: {
        name: 'A1 5\'',
        oh: 'GGAG'
      },
      a2_5: {
        name: 'A2 5\'',
        oh: 'TGAC'
      },
      a3_5: {
        name: 'A3 5\'',
        oh: 'TCCC'
      },
      b1_5: {
        name: 'B1 5\'',
        oh: 'TACT'
      },
      b2_5: {
        name: 'B2 5\'',
        oh: 'CCAT'
      },
      b3_5: {
        name: 'B3 5\'',
        oh: 'AATG'
      },
      b4_5: {
        name: 'B4 5\'',
        oh: 'AGCC',
        prev_bases: 'tc'
      },
      b5_5: {
        name: 'B5 5\'',
        oh: 'TTCG',
        prev_bases: 'tc'
      },
      b6_5: {
        name: 'B6 5\'',
        oh: 'GCTT',
        prev_bases: 'tc'
      },
      c1_5: {
        name: 'C1 5\'',
        oh: 'GGTA'
      },
      c1_3: {
        name: 'C1 3\'',
        oh: 'CGCT'
      },
    }
  },
  moclo:{
    name: 'MoClo',
    enzyme: 'bpii',
    bases_upto_snip: 'aa',
    default: {
      5: 'a',
      3: 'b'
    },
    l0_receiver: {
      name: '',
      ohs: {
        oh5: '',
        oh3: ''
      },      
    },
    domestication_enzymes:[
        'bsai',
        'bpii',
    ],
    ohs: {
      a: {
        name: 'A',
        oh: 'GGAG'
      },
      b: {
        name: 'B',
        oh: 'TACT'
      },
      c: {
        name: 'C',
        oh: 'AATG'
      },
      d: {
        name: 'D',
        oh: 'AGGT',
        prev_bases: 'tc'
      },
      e: {
        name: 'E',
        oh: 'GCTT',
        prev_bases: 'tc'
      },
      f: {
        name: 'F',
        oh: 'CGCT'
      },
    },
  }
}

export default L0DPartsStandards;
