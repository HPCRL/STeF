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
order_num=-2
order = []
times = []

pad = 6
first_time = -2
for l in f:
	if './bi' in l:
		#print('here' + l)
		tns=l.split('/')[-1].split()[0]
#        print('\n'+tns+','+str(mode),end=',')
		if len(times) > 0: # and len(order)+1 == len(times):
			first = times[0]
			times = times[1:]
			if order_num == -1:
				order = order[1:]
			sorted_times = ["N/A"]*len(order)
			for i in range(len(order)):
				x = order[i]
				sorted_times[int(x)] = times[i]
				#print(i)
				#print(times)
				#print(order)
				#print()

			for o in order:
				header += o.strip()+','

			for x in range(pad - len(order)):
				header += ','

			header += 'basic mttkrp dot time, '+first + ',ordered times ,'
			for t in sorted_times:
				header += t.strip()+','

			for x in range(pad - len(order)):
				header += ','
			if order_num == -1:
				header += 'ordered'
			print (header)
		#print()
		header=tns+','
		times = []
		order = []
		order_num=int(l.split()[-1])
		#group = (group + 1) % 4
#        if group == 0:           mode = (mode + 1) % 5
	if 'THREADS' in l:
		th = l.split('THREADS=')[-1].split()[0]
		header+= th+','
		#print('thread'+header)

#        print('th is ' +th)
#        print(l.split('THREADS='))
#        print(th,end=',')
	if 'Trying order' in l:
		order = l.split('->')[1:]
		#for x in order:            header += x.strip() + ','
	if 'its = ' in l:
		t = l.split()
		itn = int(t[2])
		if itn > it:
			it = itn

	#if 'mode' in l and 'Mode' not in l and 'thread' not in l and '==' not in l and 'reducing' not in l and 'COO' not in l:


	if 'mode' in l and 'Mode' not in l and 'thread' not in l and '==' not in l and 'reducing' not in l and 'COO' not in l:
		t = l.split()
		#print(t)
		mode = t[-2] 

		time = t[-1]
		times += [time]
		if order_num == -1:
			order += [mode]
		#print(header+str(mode)+','+time)
		#header += time + ',' 
		#print('timing '+header)
		#print()

print (header)  
print() 