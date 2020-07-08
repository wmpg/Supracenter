class Metadata:

    def __init__(self, network, code, position, name, source=None):
        
        self.network = network
        self.code = code
        self.position = position
        self.name = name
        self.source = source
        self.enabled = True

class Station:

    """
    A station object containing information and position of the station
    """

    def __init__(self, metadata, stream, edits=None, response=None):
        """
        Arguments:
            network: [String] the network of the station
            code: [String] the code of the station
            position: [position Obj] a position object of the station containing the lat/lon/elev of the station
            channel: [String] the channel of the station (BHZ, HHZ, BDF, etc.)
            name: [String] A string describing the station's name or location (for display purposes only)
            file_name: [String] Location of the .mseed file inside the working directory
        """

        self.metadata = metadata
        self.stream = stream
        self.edits = edits
        self.response = response

    def __str__(self):
        A = "Station: {:}\n".format(self.metadata.name)
        B = str(self.position) + "\n"
        try:
            C = "Ground Distance: {:}".format(self.ground_distance)
        except:
            C = ''
        return A + B + C

    def stn_distance(self, ref_pos):
        self.distance = self.metadata.position.pos_distance(ref_pos)

    def stn_ground_distance(self, ref_pos):
        self.ground_distance = self.metadata.position.ground_distance(ref_pos)