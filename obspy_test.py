
import numpy as np
import obspy

print("Reading stream")
st_real = obspy.read("F:\\Documents\\Meteor_Research\\Event\\Romania\\RO_IPH2_11.mseed")

st = obspy.read("F:\\Documents\\Meteor_Research\\Event\\Romania\\St54_4.mseed")

print("Picking first trace")
tr = st[0].copy()

# tr_real = st_real[0].copy()
# tr.station = "IPH2"
# tr.plot()
# exit()
st[0].stats.station = "IPH2"

st.write("F:\\Documents\\Meteor_Research\\Event\\Romania\\St54_4_proc.mseed")
#st[0].stats.station = "IPH2"
# tr.plot()

# print("Reading resp")
# resp = obspy.read_inventory("F:\\Documents\\Meteor_Research\\Event\\Romania\\RO_IPH2_11.xml")

# print(tr)
# # print(tr_real)
# print(resp)

# print("Remove resp")
# tr.remove_response(inventory=resp)


# tr.plot()