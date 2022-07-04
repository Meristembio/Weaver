FWD_OR_REV = (
    ('f', 'FWD'),
    ('r', 'REV'),
)

CHECK_STATES = (
    (0, 'Not required'),
    (1, 'Pending'),
    (2, 'Correct')
)

CHECK_METHODS = (
    (0, 'Digestion'),
    (1, 'Colony PCR'),
)

SEQUENCING_STATES = (
    (0, 'Not required'),
    (1, 'Required'),
    (2, 'Correct')
)

# bootstrap defaults
COLORS = (
    ('primary', 'Blue'),
    ('secondary', 'Gray'),
    ('success', 'Green'),
    ('info', 'Cyan'),
    ('warning', 'Yellow'),
    ('danger', 'Red'),
    ('light', 'White'),
    ('dark', 'Black'),
)