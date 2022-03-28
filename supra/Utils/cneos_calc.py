
import numpy as np


def cneosVectors(vx, vy, vz):


	az = np.arctan2(vx, vy)

	vh = np.sqrt(vx**2 + vy**2)

	ze = np.arctan2(vh, vz)

	az = np.degrees(az)
	ze = np.degrees(ze)

	az %= 360
	ze %= 180

	return az, ze

if __name__ == "__main__":
	print(cneosVectors(-13.2, 8.1, 1.2))


