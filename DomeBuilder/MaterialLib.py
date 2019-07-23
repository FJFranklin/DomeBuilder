class MaterialLib(object):

    __steel = None

    def __init__(self):
        # TODO: Add materials from file
        None

    @staticmethod
    def steel():
        if MaterialLib.__steel is None:
            MaterialLib.__steel = {}
            MaterialLib.__steel['Young']   = 209E9
            MaterialLib.__steel['Shear']   = 79.3E9
            MaterialLib.__steel['density'] = 7850
            MaterialLib.__steel['yield']   = 200E6
        return MaterialLib.__steel
