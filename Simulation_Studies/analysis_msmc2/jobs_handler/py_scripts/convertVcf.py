import sys


def main():

	inFile = sys.argv[1]

	ifs = open (inFile)

	doubleHit = 0
	previousPos = 0
	for line in ifs:
		if (line.startswith("#")):
			print (line[:-1])
		else:
			fields = line.split()
			# print (fields)
			pos = int(fields[1])
			fields[3] = 'A'
			fields[4] = 'T'
			if (pos != previousPos):
				print ("\t".join(fields))
			else:
				doubleHit += 1
			previousPos = pos
	# print (doubleHit)


if __name__ == "__main__":
	main()