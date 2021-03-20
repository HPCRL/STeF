import sys


def generate_code(codes):
	tab = 0;
	for snippet in codes:
		for line in snippet:
			if '}' in line:
				tab -= 1
			for t in range(tab):
				print("\t",end='')
			print(line)
			if '{' in line:
				tab += 1


def generate_code5(codes,n):
	tab = 0
	start = codes[0]
	repeat_1 = codes[1]
	repeat_2 = codes[3]
	inside = codes[2]
	end = codes[4]

	#write the start
	for line in start:
		if '}' in line:
			tab -= 1
		for t in range(tab):
			print("\t",end='')
		print(line.replace('NNN',str(n)))
		if '{' in line:
			tab += 1


	# write the first repeat loop
	for x in range(n-2):
		for line in repeat_1:
			if '}' in line:
				tab -= 1
			for t in range(tab):
				print("\t",end='')
			print(line.replace('XXX',str(x)).replace('YYY',str(x+1)))
			if '{' in line:
				tab += 1

	#write inside loop
	for line in inside:
		if '}' in line:
			tab -= 1
		for t in range(tab):
			print("\t",end='')
		print(line.replace('XXX',str(n-2)).replace('YYY',str(n-1)))
		if '{' in line:
			tab += 1


	# write the second repeat loop
	for x in reversed(range(n-2)):
		#print ('x is '+str(x))
		for line in repeat_2:
			if '}' in line:
				tab -= 1
			for t in range(tab):
				print("\t",end='')
			print(line.replace('XXX',str(x)).replace('YYY',str(x+1)))
			if '{' in line:
				tab += 1


	#write the end
	for line in end:
		if '}' in line:
			tab -= 1
		for t in range(tab):
			print("\t",end='')
		print(line.replace('NNN',str(n)))
		if '{' in line:
			tab += 1



def readfile(file):
	f = open(file).readlines()
	res = []
	for l in f:
		res += [l.replace('\n','')]
		#print(l[:-1])

	#l = l.replace('\\','\\\\')

	return res

first_start=['int mttkrp_hardwired_first_NNN(csf* t, int mode, int r, matrix** mats, int profile )',
	'{',
	"int nmode = t->nmode;",
	"int num_th = 1;",
	"int partial_results_size = nmode*r+PAD;",
	"#ifdef OMP",
	"num_th = omp_get_max_threads();",
	"#endif",

	'printf("num ths %d\\n", num_th);',

	'TYPE* partial_results_all = (TYPE*) malloc(num_th*partial_results_size*sizeof(TYPE));',

	'for(int i = 0 ; i < num_th*partial_results_size ; i++)',
		'partial_results_all[i] = 0;',

	'set_matrix(*mats[0],0);',

	'#ifdef OMP',
	'#pragma omp parallel ',
	'#endif',
	'{',
		'int th = 0;',
		'#ifdef OMP',
		'th = omp_get_thread_num();',
		'#endif',
		'auto time_start = std::chrono::high_resolution_clock::now();',
		'TYPE* partial_results = partial_results_all + th*partial_results_size;',
		'#ifdef OMP',
		'#pragma omp for schedule(dynamic,1)',
		'#endif',
		'for(idx_t i0 = 0 ; i0< t->fiber_count[0]; i0++)',
		"{"]

first_repeat_1=['for(idx_t iYYY = t->ptr[XXX][iXXX] ; iYYY < t->ptr[XXX][iXXX+1]; iYYY++)',
			'{'	]

first_inside=['for(idx_t iYYY = t->ptr[1][i1] ; iYYY< t->ptr[XXX][iXXX+1]; iYYY++)',
				'{',
					'TYPE* pr = partial_results + XXX * r;',
					'TYPE tval = t->val[iYYY];',
					'TYPE* matval = (mats[YYY]->val) + ((mats[YYY]) -> dim2) * t->ind[YYY][iYYY];',
					'#pragma omp simd',
					'for(int y=0 ; y<r ; y++)',
					'{',
						'pr[y] += tval * matval[y];	// TTM step			',
				'}' , '}']

first_repeat_2=['TYPE* matval = (mats[YYY]->val) + ((mats[YYY]) -> dim2) * t->ind[YYY][iYYY];',
				'TYPE* intval = t->intval[YYY] + iYYY*r*YYY;',
				'#pragma omp simd',
				'for(int y=0 ; y<r ; y++)',
				'{',
					#'//printf("1st level loop %lf\n",partial_results[r+y]);',
					'partial_results[y] += partial_results[r+y] * matval[y]; // TTV',
					'intval[y] = partial_results[YYY*r + y];',
					
				'}',
				'#pragma omp simd',
				'for(int y=0; y<r; y++)',
				'{',
					'partial_results[YYY*r+y] = 0;',
				'}', '}'	]

first_end=['// write to output matrix',
			'TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0];',

			'#pragma omp simd',
			'for(int y=0 ; y<r ; y++)',
			'{',
				
				'matval[y] = partial_results[y];',
				'partial_results[y] = 0;',
			'}',

		'}',

		'auto time_end = std::chrono::high_resolution_clock::now();',
		'std::chrono::duration<double> time_diff = time_end-time_start;',
		
		'printf("Hardwired time for mode %d thread %d %lf \\n",t->modeid[mode],th,time_diff.count());',

	"}",
	'rem(partial_results_all);	',
	'return 0;', '}']



'''
start = open(sys.argv[1]).readlines()
r1 = open(sys.argv[2]).readlines()
ins = open(sys.argv[3]).readlines()
r2 = open(sys.argv[4]).readlines()
end = open(sys.argv[5]).readlines()
'''
code_order = [first_start,first_repeat_1,first_inside,first_repeat_2,first_end]
#generate_code(code_order)

#generate_code5(code_order,5)
#for x in range(10):
#	print('###################################	')

#co = [start,r1,ins,r2,end]
co = [readfile(x) for x in sys.argv[1:6]]
dim=5
if len(sys.argv)>6:
	dim = int(sys.argv[6])
generate_code5(co,dim)