def str_time(seconds):
	res = ""
	if seconds < 60:
		res = str(seconds) + " s"
	elif seconds < 3600:
		res = str(round(seconds / 60.0, 2)) + " min"
	elif seconds < 3600 * 24:
		res = str(round(seconds / 3600.0, 2)) + " h"
	elif seconds < 3600 * 24 * 7:
		res = str(round(seconds / (3600 * 24.0), 2)) + " d"
	else:
		res = str(seconds) + " s"

	return res