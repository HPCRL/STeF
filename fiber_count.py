import sys

f = open(sys.argv[1]).readlines()

th=-1
tns="tensor"
it=-1
sumind=0


accept=['L1','L2','DRAM','FP_ARITH']
accept+=['BR','STALL','Branch','Stall']

header=''
mode = 0
group = 0
event=False
first_time=-1
order_num=-2
order = []
times = []
modes = []
fiber = []
pad=6

for l in f:
	if './bi' in l:
		#print('here' + l)
		tns=l.split('/')[-1].split()[0]
#        print('\n'+tns+','+str(mode),end=',')
		#if order_num == -1:

		order = order[1:]

			#print(order)
			#print(times)
		header += ' mode order,'
		for o in order:
			header += o+','
		for x in range(pad - len(order)):
			header += ','
		header += ' mode lengths,'
		for t in modes:
			header += t+','
		for x in range(pad - len(order)):
			header += ','
		header += ' fiber lengths,'
		for t in fiber:
			header += t+','
		for x in range(pad - len(order)):
			header += ','
		
		if order_num == -1:
			header += 'ordered'
		if len(order) > 0:
			print (header)
		#print()
		header=tns+','
		times = []
		order = []
		modes = []
		fiber = []
		order_num=int(l.split()[-1])
		#if order_num == -1:			print(header)
		#group = (group + 1) % 4
#        if group == 0:           mode = (mode + 1) % 5
	if 'THREADS' in l:
		th = l.split('THREADS=')[-1].split()[0]
		header+= th+','

		#print('thread'+header)
		#if order_num == -1:			print('th is ' +th)
#        print(l.split('THREADS='))
#        print(th,end=',')

	if 'its = ' in l:

		t = l.split()
		itn = int(t[2])
		if itn > it:
			it = itn
	'''
	if 'no fusion' in l and False:
		t = l.split()
		time = t[-1]
		if order_num == -1:
			times += [time]
		else:
			header += time + ',' 
	'''

	if "Mode" in l:
		t = l.split()
		#print(l)
		if 'fiber' in l:
			fiber += [t[-1]]
		if 'mode' in l:
			modes += [t[-1]]
	
	if 'mode' in l and 'Mode' not in l and 'thread' not in l and '==' not in l and 'reducing' not in l and 'COO' not in l:
		t = l.split()
		#print(t)
		mode = t[-2] 

		time = t[-1]

		order += [mode]


#print (header)	
print() 

