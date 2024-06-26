import re


# the regex to find the lines we want from the report
prog = re.compile(r'[\d-]+\s+[a-zA-Z():\/]*\s\d+\s+[a-zA-Z():\/\s.0-9 -]*')

units = ['Date', 'Clearance', 'Normalized Clearance', 'GFR']

def read_report(filename):
	with open (filename) as file:
		lines = file.readlines()
	return lines


def get_data(lines):
	for line in lines:
		data = prog.match(line)
		if data != None:
			split_line = re.split(r'\s\s\s\s', data.group(0))
			split_line = [l.strip() for l in split_line]


			split_line[1] = re.findall(r'[0-9]+', split_line[1])[0]
			split_line[2] = re.findall(r'\d+\s\([\d\s-]+\)', split_line[2])[0]
			print(dict(zip(units,split_line)))


if __name__ == '__main__':
	text = read_report('/run/user/1000/gvfs/smb-share:server=rghheapac001.regionh.top.local,share=data/TEMP/GFR_test/pdftest/REPORTS/0000B695/AA1358F5/AA96E1D7.RTF' 
)

	get_data(text)



