# add your own OH standard
standards = {
    'loop': {
        'odd': (
            ('A', 'GGAG'),
            ('B', 'TACT'),
            ('C', 'AATG'),
            ('X', 'TGGA'),
            ('D', 'AGGT'),
            ('E', 'GCTT'),
            ('F', 'CGCT'),
            ('None', ''),
        ),
        'odd_custom': (
            ('A', 'GGAG'),
            ('B', 'TACT'),
            ('C', 'AATG'),
            ('X', 'TGGA'),
            ('D', 'AGGT'),
            ('E', 'GCTT'),
            ('F', 'CGCT'),
            ('None', ''),
        ),
        'even': (
            ('ALPHA', 'ATG'),
            ('BETA', 'GCA'),
            ('GAMMA', 'TAC'),
            ('EPSILON', 'CAG'),
            ('OMEGA', 'GGT'),
        )
    }
}

# Modify to your needs
CURRENT_ASSEMBLY_STANDARD = standards['loop']
