"""Script to download atmospheric data file from NASA Earthdata MERRA-2"""

def fetchMERRA(setup):
    """ Python code taken from: https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
        to fetch data from NASA Earthdata MERRA-2. Saves the weather file in file_name. Only downloads the 
        required part of the data, usually only a few megabytes.

    Arguments:
        file_name: [string] location to save the file
        year, month, day, hour: [int] time to take weather data from
        search_area: [list] boundary to take weather data from [lat min, lat max, lon min, lon max]
        username: [string] NASA Earthdata username
        password: [string] NASA Earthdata password
        
    """

    # search a little around the boundary for weather profiles for positions close to the edge
    lat_min = setup.search_area[0] - 0.25
    lat_max = setup.search_area[1] + 0.25 
    lon_min = setup.search_area[2] - 0.25
    lon_max = setup.search_area[3] + 0.25

    #convert to strings
    y = str(setup.ref_time.year)
    m = str(setup.ref_time.month).zfill(2)
    d = str(setup.ref_time.day).zfill(2)
    h = str(setup.ref_time.hour).zfill(2)
    lat_m = str(lat_min)
    lat_M = str(lat_max)
    lon_m = str(lon_min)
    lon_M = str(lon_max)

    username = setup.username
    password = setup.password

    url = '''http://goldsmr5.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FMERRA2%2FM2I3NVASM.5.12.4%2F'''+y+'''%2F'''+m+'''%2FMERRA2_400.inst3_3d_asm_Nv.'''+y+m+d+'''.nc4&FORMAT=aGRmLw&BBOX='''+lat_m+'''%2C'''+lon_m+'''%2C'''+lat_M+'''%2C'''+lon_M+'''&TIME='''+y+'''-'''+m+'''-'''+d+'''T'''+h+'''%3A00%3A00%2F'''+y+'''-'''+m+'''-'''+d+'''T'''+h+'''%3A00%3A00&LABEL=MERRA2_400.inst3_3d_asm_Nv.'''+y+m+d+'''.SUB.hdf&SHORTNAME=M2I3NVASM&SERVICE=SUBSET_MERRA2&VERSION=1.02&DATASET_VERSION=5.12.4&VARIABLES=H%2CT%2CU%2CV'''
     
    import requests # get the requsts library from https://github.com/requests/requests
     
    # overriding requests.Session.rebuild_auth to mantain headers when redirected
    class SessionWithHeaderRedirection(requests.Session):
     
        AUTH_HOST = 'urs.earthdata.nasa.gov'
     
        def __init__(self, username, password):
     
            super().__init__()
     
            self.auth = (username, password)
     
       # Overrides from the library to keep headers when redirected to or from
       # the NASA auth host.
        def rebuild_auth(self, prepared_request, response):
     
            headers = prepared_request.headers
     
            url = prepared_request.url

            if 'Authorization' in headers:
     
                original_parsed = requests.utils.urlparse(response.request.url)
     
                redirect_parsed = requests.utils.urlparse(url)
     
                if (original_parsed.hostname != redirect_parsed.hostname) and \
                        redirect_parsed.hostname != self.AUTH_HOST and \
                        original_parsed.hostname != self.AUTH_HOST:
     
                    del headers['Authorization']

            return
     
    # create session with the user credentials that will be used to authenticate access to the data
    session = SessionWithHeaderRedirection(username, password)
     
    # extract the filename from the url to be used when saving the file
    # filename = url[url.rfind('/')+1:]  

    try:
     
        # submit the request using the session  
        response = session.get(url, stream=True)
     
        # raise an exception in case of http errors
        response.raise_for_status()  
     
        # save the file
        with open(file_name, 'wb') as fd:
     
            for chunk in response.iter_content(chunk_size=1024*1024):
     
                fd.write(chunk)
     
      
     
    except requests.exceptions.HTTPError as e:
     
        # handle any errors here
        print(e)