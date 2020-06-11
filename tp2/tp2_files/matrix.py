def matrix_print(a, h, w) :
	print ("matrix is a %dx%d" % (h, w))
	i = 0
	while (i < h) :
		j = 0
		while (j < w) :
			print ("%d\t" % a[i * w + j]),
			j += 1
		print ""
		i += 1

def matrix_print_sub(a, h, w, k, l, bh, bw) :
	print ("sub matrix is a %dx%d, line %d column %d" % (bh, bw, k, l))
	i = 0
	while (i < bh) :
		j = 0
		while (j < bw) :
			print ("%d\t" % a[(i + k) * w + (j + l)]),
			j += 1
		print ("")
		i += 1

def main() :
	# array (representing a 2d matrix)
	m = [
		11, 12, 13, 13, 14, 15, 16, 17,
		21, 22, 23, 24, 25, 26, 27, 28,
		31, 32, 33, 34, 35, 36, 37, 38,
		41, 42, 43, 44, 45, 46, 47, 48,
		51, 52, 53, 54, 55, 56, 57, 58,
		61, 62, 63, 64, 65, 66, 67, 68,
		71, 72, 73, 74, 75, 76, 77, 78,
		81, 82, 83, 84, 85, 86, 87, 88,
		91, 92, 93, 94, 95, 96, 97, 98
	]
	# print the matrix
	matrix_print(m, 9, 8)
	# matrix, height, width, block line, block column, block height, block width
	matrix_print_sub(m, 9, 8, 1, 2, 3, 2)

if __name__ == "__main__" : main()
