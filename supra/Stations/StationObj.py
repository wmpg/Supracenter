class Metadata:

    def __init__(self, network, code, position, name, source=None):
        
        self.network = network
        self.code = code
        self.position = position
        self.name = name
        self.source = source
        self.enabled = True

class Polarization:

    def __init__(self):

        self.azimuth = []
        self.azimuth_error = []
        self.time = []

    def __str__(self):
        return "Polarization object with azimuth = {:} +/- {:}".format(self.azimuth, self.azimuth_error)


class AnnotationList:

    def __init__(self):

        self.annotation_list = []

    def __str__(self):

        A = ""

        if len(self.annotation_list) > 0:
            for item in self.annotation_list:
                A += "{:}\n".format(item.title)
        else:
            A = "No annotations"

        return A

    def add(self, an):

        self.annotation_list.append(an)

    def overwrite(self, an):

        new_list = []

        for annote in self.annotation_list:
            if annote.title != an.title:
                new_list.append(annote)

        new_list.append(an)

        self.annotation_list = new_list

    def remove(self, an):
        
        new_list = []

        for annote in self.annotation_list:
            if annote.title != an.title:
                new_list.append(annote)

        self.annotation_list = new_list


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
        self.polarization = Polarization()
        self.annotation = AnnotationList()
        self.bandpass = None

    def __str__(self):
        A = "Station: {:}\n".format(self.metadata.name)

        B = str(self.metadata.position) + "\n"

        try:
            C = "Ground Distance: {:} \n".format(self.ground_distance)
        except:
            C = ''

        if self.response is not None:
            D = "Response Attatched"
        else:
            D = "No Response Found"


        return A + B + C + D

    def hasResponse(self):

        if not hasattr(self, "response"):
            return False

        if self.response is None:
            return False

        return True

    def hasSeismic(self):

        st = self.stream
        a = st.select(component="Z")
        if len(a) > 0:
            return True
        else:
            return False

    def hasInfrasound(self):
        st = self.stream
        a = st.select(component="F")
        if len(a) > 0:
            return True
        else:
            return False

    def stn_distance(self, ref_pos):
        self.distance = self.metadata.position.pos_distance(ref_pos)
        return self.distance

    def stn_ground_distance(self, ref_pos):
        self.ground_distance = self.metadata.position.ground_distance(ref_pos)