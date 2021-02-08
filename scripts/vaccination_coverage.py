import numpy as np
from collections import defaultdict

def read_country_codes(fname):
	cc = {}
	with open(fname) as fh:
		for line in fh:
			if line[0]=='#':
				continue
			try:
				c, country = line.strip().split('\t')
				cc[c]=country
			except:
				print(line)
	return cc

def read_coverage(fname):
	vaccov = {}
	with open(fname) as fh:
		for line in fh:
			if line[0]=='#':
				continue
			try:
				c, cov = line.strip().split('\t')
				vaccov[c]=float(cov)
			except:
				print(line)
	return vaccov


def read_OECD_coverage(fname, average=None, country_codes=None):
	vaccov = defaultdict(list)
	if country_codes is None:
		country_codes = {}
	with open(fname, 'r') as fh:
		header = fh.readline()
		for line in fh:
			entries = line.strip().split(',')
			vaccov[country_codes.get(entries[0],entries[0])].append((int(entries[5]), float(entries[6])))

	for c, d in vaccov.items():
		if average:
			vaccov[c] = np.mean(np.array(sorted(d, key=lambda x:x[0]), dtype=float)[-average:,1])
		else:
			vaccov[c] = np.array(sorted(d, key=lambda x:x[0]), dtype=float)

	return vaccov


def read_all_vaccination_data():
	country_codes = read_country_codes('source-data/country_codes.tsv')
	vaccov_OECD = read_OECD_coverage('source-data/2018_OECD_flu_vaccination_coverage.csv', average=3, country_codes=country_codes)
	vaccov_europe = read_coverage('source-data/2018_Europe_flu_vaccination_coverage.tsv')
	vaccov_SA = read_coverage('source-data/2018_South_America_flu_vaccination_coverage.tsv')

	vaccov = defaultdict(list)
	for d in [vaccov_SA, vaccov_europe, vaccov_OECD]:
		for c,v in d.items():
			vaccov[c].append(v)
	for c,v in vaccov.items():
		vaccov[c] = np.median(v)

	return vaccov

if __name__ == '__main__':
	vaccov = read_all_vaccination_data()
