import sys

TO_EXCLUDE = [
	"PROBLEM CLASS", 
	"N. OF ITEMS", 
	"RELATIVE AND ABSOLUTE N. OF INSTANCE",
	"HBIN,WBIN",
	"H(I),W(I),I=1,...,N" 
]

def read_elem(filename):
	data = []
	with open(filename) as f:
		for line in f:
			for word in TO_EXCLUDE:
				if word in line:
					line = line.replace(word, "")
			line = line.strip()
			line = " ".join(line.split())
			data.append(line)
	return data

def read_input_b2p(filename):
	raw = read_elem(filename)
	count = 0
	n = 0
	data = []
	bins = []

	for el in raw:
		if not el.strip():
			continue
		else:
			if count == 0:
				problem_class = int(el)
				count += 1 
			elif count == 1:
				n_elements = int(el)
				count += 1 
			elif count == 2:
				pos_relative, pos_absolute = [int(i) for i in el.split(" ")]
				count += 1 
			elif count == 3:
				hbin, wbin =  [int(i) for i in el.split(" ")]
				count += 1 

			else:
				if n < n_elements-1:
					bins.append([int(i) for i in el.split(" ")])
					n += 1
				else:
					bins.append([int(i) for i in el.split(" ")])
					data.append((bins, hbin, wbin, problem_class))
					count =  0
					bins = []
					n = 0

	return data

if __name__ == '__main__':
	data = read_input_b2p("Instancias/BW_2bp/Class_01.2bp")
	assert len(data) == 50