const AssemblyStandards = {
    "ytk": {
        name:"YTK",
        ohs: {
                '81': {
                    'name': '8-1',
                    'oh': 'CCCT'
                },
                '12': {
                    'name': '1-2',
                    'oh': 'AACG'
                },
                '23': {
                    'name': '2-3',
                    'oh': 'TATG'
                },
                '3a3b': {
                    'name': '3a-3b',
                    'oh': 'TTCT',
                    'prev_bases': 'gg'
                },
                '34': {
                    'name': '3-4',
                    'oh': 'ATCC'
                },
                '4a4b': {
                    'name': '4a-4b',
                    'oh': 'TGGC'
                },
                '45': {
                    'name': '4-5',
                    'oh': 'GCTG'
                },
                '56': {
                    'name': '5-6',
                    'oh': 'TACA'
                },
                '67': {
                    'name': '6-7',
                    'oh': 'GAGT'
                },
                '78': {
                    'name': '7-8',
                    'oh': 'CCGA'
                },
                '8a8b': {
                    'name': '8a-8b',
                    'oh': 'CAAT'
                },
                'S': {
                    'name': 'S',
                    'oh': 'CTGA'
                },
                '1': {
                    'name': '1',
                    'oh': 'CCAA'
                },
                '2': {
                    'name': '2',
                    'oh': 'GATG'
                },
                '3': {
                    'name': '3',
                    'oh': 'GTTC'
                },
                '4': {
                    'name': '4',
                    'oh': 'GGTA'
                },
                '5': {
                    'name': '5',
                    'oh': 'AAGT'
                },
                'E': {
                    'name': 'E',
                    'oh': 'AGCA'
                },
        },
        enzymes: [
            {
                value: 'BsaI',
                name: 'L1 (BsaI)'
            },
            {
                value: 'BsmBI',
                name: 'L2 (BsmBI)'
            },
        ]
    },
    "loop": {
        name:"Loop Assembly",
        ohs: {
            a: {
                'name': 'A',
                'oh': 'GGAG'
            },
            b: {
                'name': 'B',
                'oh': 'TACT'
            },
            c: {
                'name': 'C',
                'oh': 'AATG'
            },
            d: {
                'name': 'D',
                'oh': 'AGGT',
            },
            e: {
                'name': 'E',
                'oh': 'GCTT'
            },
            f: {
                'name': 'F',
                'oh': 'CGCT'
            },
            alpha: {
                'name': 'Alpha',
                'oh': 'ATG'
            },
            beta: {
                'name': 'Beta',
                'oh': 'GCA'
            },
            Gamma: {
                'name': 'Gamma',
                'oh': 'TAC'
            },
            Epsilon: {
                'name': 'Epsilon',
                'oh': 'CAG'
            },
            Omega: {
                'name': 'Omega',
                'oh': 'GGT'
            },
        },
        enzymes: [
            {
                value: 'BsaI',
                name: 'ODD (BsaI / Eco31I)'
            },
            {
                value: 'SapI',
                name: 'EVEN (SapI / LguI)'
            },
        ]
    },
    "moclo": {
        name:"MoClo",
        ohs: {
            p5: {
                'name': 'P.5\'',
                'oh': 'GGAG'
            },
            p3u5: {
                'name': 'P.3\'-U.5\'',
                'oh': 'TACT'
            },
            u3sp5: {
                'name': 'U.3\'-SP.5\'',
                'oh': 'AATG'
            },
            sp3cds5: {
                'name': 'SP.3\'-CDS.5\'',
                'oh': 'AGGT'
            },
            cds3t5: {
                'name': 'CDS.3\'-T.5\'',
                'oh': 'GCTT'
            },
            t3: {
                'name': 'T.3\'',
                'oh': 'CGCT'
            },
            l115l173: {
                'name': 'L1-1.5\'-L1-7.3\'',
                'oh': 'TGCC'
            },
            l125l113: {
                'name': 'L1-2.5\'-L1-1.3\'',
                'oh': 'GCAA'
            },
            l135l123: {
                'name': 'L1-3.5\'-L1-2.3\'',
                'oh': 'ACTA'
            },
            l145l133: {
                'name': 'L1-4.5\'-L1-3.3\'',
                'oh': 'TTAC'
            },
            l155l143: {
                'name': 'L1-5.5\'-L1-4.3\'',
                'oh': 'CAGA'
            },
            l165l153: {
                'name': 'L1-6.5\'-L1-5.3\'',
                'oh': 'TGTG'
            },
            l175l163: {
                'name': 'L1-7.5\'-L1-6.3\'',
                'oh': 'GAGC'
            },
        },
        enzymes: [
            {
                value: 'BsaI',
                name: 'L1 (BsaI / Eco31I)'
            },
            {
                value: 'BbsI',
                name: 'L2 (BpiI / BbsI)'
            },
        ]
    },
    "generic": {
        name:"Generic",
        ohs: {},
        enzymes: [        
            {
                value: 'AarI',
                name: 'AarI / PaqCI'
            },
            {
                value: 'BbsI',
                name: 'BbsI / BpiI'
            },
            {
                value: 'BsaI',
                name: 'BsaI / Eco31I'
            },
            {
                value: 'BsmBI',
                name: 'BsmBI / Esp3I'
            },
            {
                value: 'SapI',
                name: 'LguI / SapI'
            },
        ]
    }
}

export default AssemblyStandards;