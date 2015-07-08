#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
Mathieu Daëron, May 2012 (mathieu@daeron.fr)
Last tested using Python 2.7.1
This code is licensed under a Creative Commons
Attribution-NonCommercial-ShareAlike 3.0 Unported License.
[http://creativecommons.org/licenses/by-nc-sa/3.0/]
'''

import matplotlib as mpl
mpl.rcParams['backend'] = 'Agg'
mpl.rcParams['savefig.dpi'] = '300'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial']
mpl.rcParams['font.size'] = '10'
mpl.rcParams['axes.labelsize'] = '10'
mpl.rcParams['pdf.fonttype'] = '42'

import numpy, pylab
from pylab import linalg

# Reduced data from the mass spectrometer;
# lines starting with 'e' correspond to CO2 equilibrated to a known temperature (Teq);
# lines starting with 'u' correspond to unknowns;
# 'sD47' is the (internal) standard error for each measurement;
# 'id' stands for in-lab tracking numbers.

ms_data = """type	Teq	d47	D47	sD47	id	label
e	1000	.1	-.630	.008	1	
e	1000	.2	-.621	.010	2	
e	1000	.3	-.642	.011	3	
e	1000	.4	-.635	.009	4	
e	1000	.5	-.629	.010	5	
e	1000	19.1	-.430	.008	6	
e	1000	19.2	-.421	.010	7	
e	1000	19.3	-.442	.011	8	
e	1000	19.4	-.435	.009	9	
e	1000	19.5	-.429	.010	10	
e	20	.1	.220	.008	11	
e	20	.2	.221	.010	12	
e	20	.3	.242	.011	13	
e	20	.4	.235	.009	14	
e	20	.5	.229	.010	15	
e	20	19.1	.520	.008	16	
e	20	19.2	.521	.010	17	
e	20	19.3	.542	.011	18	
e	20	19.4	.535	.009	19	
e	20	19.5	.529	.010	20	
u		15.2	.095	.010	24	Sample A
u		13.0	-.054	.006	25	Sample B
u		9.50	-.225	.012	26	Sample C"""

# Theoretical equilibrium values for Δ47 in CO2 gas according to
# Wang et al. (2004) [http://dx.doi.org/10.1016/j.gca.2004.05.039],
# reported in the supplementary data of
# Dennis et al. (2011) [http://dx.doi.org/10.1016/j.gca.2011.09.025]:

eq_val = '''-83	1.8954
-73	1.7530
-63	1.6261
-53	1.5126
-43	1.4104
-33	1.3182
-23	1.2345
-13	1.1584
-3	1.0888
7	1.0251
17	0.9665
27	0.9125
37	0.8626
47	0.8164
57	0.7734
67	0.7334
87	0.6612
97	0.6286
107	0.5980
117	0.5693
127	0.5423
137	0.5169
147	0.4930
157	0.4704
167	0.4491
177	0.4289
187	0.4098
197	0.3918
207	0.3747
217	0.3585
227	0.3431
237	0.3285
247	0.3147
257	0.3015
267	0.2890
277	0.2771
287	0.2657
297	0.2550
307	0.2447
317	0.2349
327	0.2256
337	0.2167
347	0.2083
357	0.2002
367	0.1925
377	0.1851
387	0.1781
397	0.1714
407	0.1650
417	0.1589
427	0.1530
437	0.1474
447	0.1421
457	0.1370
467	0.1321
477	0.1274
487	0.1229
497	0.1186
507	0.1145
517	0.1105
527	0.1068
537	0.1031
547	0.0997
557	0.0963
567	0.0931
577	0.0901
587	0.0871
597	0.0843
607	0.0816
617	0.0790
627	0.0765
637	0.0741
647	0.0718
657	0.0695
667	0.0674
677	0.0654
687	0.0634
697	0.0615
707	0.0597
717	0.0579
727	0.0562
737	0.0546
747	0.0530
757	0.0515
767	0.0500
777	0.0486
787	0.0472
797	0.0459
807	0.0447
817	0.0435
827	0.0423
837	0.0411
847	0.0400
857	0.0390
867	0.0380
877	0.0370
887	0.0360
897	0.0351
907	0.0342
917	0.0333
927	0.0325
937	0.0317
947	0.0309
957	0.0302
967	0.0294
977	0.0287
987	0.0281
997	0.0274
1007	0.0268
1017	0.0261
1027	0.0255
1037	0.0249
1047	0.0244
1057	0.0238
1067	0.0233
1077	0.0228
1087	0.0223
1097	0.0218'''

eq_val = numpy.array([l.split('\t') for l in eq_val.split('\n')], dtype='float')
T = eq_val[:,0] - 0.15
D47 = eq_val[:,1]

# Quadratic polynomial fit as used by
# Dennis et al. (2011) [http://dx.doi.org/10.1016/j.gca.2011.09.025]:
eq_fit = numpy.polyfit( 1/(T+273.15), D47, 4 )
EqD47 = lambda T : \
	eq_fit[0] * (T+273.15)**-4 + \
	eq_fit[1] * (T+273.15)**-3 + \
	eq_fit[2] * (T+273.15)**-2 + \
	eq_fit[3] * (T+273.15)**-1 + \
	eq_fit[4]




# Read the MS data:
E = []
S = []

for l in ms_data.split('\n')[1:] :
	l = l.split('\t')
 	if l[0] == 'e' :
 		Teq = float(l[1])
 		if Teq not in [e['Teq'] for e in E] :
 			E.append({ 'Teq':Teq, 'D47eq':EqD47(Teq), 'data':[] })
		for e in E :
			if e['Teq'] == Teq :
				e['data'].append({
					'd47':float(l[2]),
					'D47':float(l[3]),
					'sD47':float(l[4]),
					'id':l[5],
					'label':l[6]
					})
	if l[0] == 'u' :
		S.append({
			'd47':float(l[2]),
			'D47':float(l[3]),
			'sD47':float(l[4]),
			'id':l[5],
			'label':l[6]
			})

# Define function to fit the 'E' data
# to a plane in (d47, rawD47, theoreticalD47) space:

def fit_allgasdata( E, conservative=True ) :
	A = [] # design matrix
	Y = [] # target values for the fit
	for e in E :
		for d in e['data'] :
			A.append( [ e['D47eq']/d['sD47'], d['d47']/d['sD47'], 1./d['sD47'] ] )
			Y.append( d['D47'] / d['sD47'] )
	A,Y = numpy.array(A), numpy.array(Y)
	f = linalg.lstsq(A,Y.T)[0] # best-fit parameters
	CM = linalg.inv(linalg.dot(A.T,A)) # covariance matrix of fit parameters
	if conservative :
		# Scale up uncertainties in the fit parameters if the goodnes-of-fit is worse than average.
		# To some extent, this helps account for external errors in the gas line data.
		chi2 = sum( ( Y - linalg.dot( A, f ) )**2 )
		nf = sum([len(e['data']) for e in E]) -3
		if chi2 > nf :
			CM = CM * chi2 / nf
	return f,CM

(a,b,c),CM = fit_allgasdata(E)

print 'id	D47	sD47	label'
for s in S :

	x = s['D47']
	y = s['d47']
	sx = s['sD47']
	sy = 0
	z = x / a - b/a * y - c/a
	s['trueD47'] = z	
	
	dzdx = a ** -1
	dzdy = -b / a
	dzda = -(x-b*y-c) * a ** -2
	dzdb = -y / a
	dzdc = -a ** -1

	v = numpy.array([ dzdx, dzdy, dzda, dzdb, dzdc ])
	C = numpy.zeros((5,5))
	C[0,0] = sx ** 2
	C[1,1] = sy ** 2
	C[2:,2:] = CM
	sz = (linalg.dot( v, linalg.dot( C, v.T ) ))**.5
	s['strueD47'] = sz
 
	# Although the above formula holds true for linear combinations of normal variables
	# (Tellinghuisen, 2001) [http://dx.doi.org/10.1021/jp003484u],
	# Our Monte-Carlo simulations showed that distributions of D47 values
	# estimated using this code are undistinguishable from normal distributions,
	# so that 95% confidence limits correspond to +/- 2*sD47corrected.
	# These limits account for internal and external uncertainties on gas line data,
	# and for internal errors in sample measurements, but not for external uncertainties in the latter.
	# It is thus always recommended to check external reproducibility for
	# duplicate measurements of homogeneous samples.
	
	print '%s	%.4f	%.4f	%s' %( s['id'], s['trueD47'], s['strueD47'], s['label'] )

# xmin = numpy.floor(
# 	min(
# 		min( [s['d47'] for s in S] ),
# 		min( [min([d['d47'] for d in e['data']]) for e in E] ) )
# 		) - 1.
# xmax = numpy.ceil(
# 	max(
# 		max( [s['d47'] for s in S] ),
# 		max( [max([d['d47'] for d in e['data']]) for e in E] ) )
# 		) + 1.
# xi = numpy.array([ xmin, xmax ])
# 
# for k in range(len(E)) :
# 	yi = f[0]*xi+f[k+1]
# 	pylab.plot(xi,yi,'r-')
# 	pylab.text(
# 		xi[-1],
# 		yi[-1],
# 		u' %g °C' %(E[k]['Teq']),
# 		horizontalalignment='left',
# 		verticalalignment='center',
# 		color='r',
# 		alpha=1
# 		)
# 
# Dmin = min( E[1]['D47eq'], E[0]['D47eq'] )
# Dmax = max( E[1]['D47eq'], E[0]['D47eq'] )
# Drange = numpy.r_[numpy.ceil(Dmin*10)/10:numpy.ceil(Dmax*10)/10:.1]
# 
# for D in Drange[::2] :
# 	yi = f[0] * xi + f[1] + (f[2]-f[1])*(D-E[0]['D47eq'])/(E[1]['D47eq']-E[0]['D47eq'])
# 	pylab.plot(xi,yi,'r-',alpha=.25)
# 
# for D in Drange[1:-1:2] :
# 	yi = f[0] * xi + f[1] + (f[2]-f[1])*(D-E[0]['D47eq'])/(E[1]['D47eq']-E[0]['D47eq'])
# 	pylab.plot(xi,yi,'r-',alpha=.25)
# 	pylab.text(
# 		xi[-1],
# 		yi[-1],
# 		u' %.1f ‰' %(D),
# 		horizontalalignment='left',
# 		verticalalignment='center',
# 		color='r',
# 		alpha=.5
# 		)
# 
# for e in E :
# 	for d in e['data'] :
# 		pylab.errorbar(d['d47'],d['D47'],2*d['sD47'],None,'wo',mew=1,mec='r',ecolor='r',ms=5,capsize=0)
# for s in S :
# 	pylab.errorbar(s['d47'],s['D47'],2*s['sD47'],None,'wo',mew=1,mec='k',ecolor='k',ms=5,capsize=0)
# 	pylab.text(
# 		s['d47'],
# 		s['D47'],
# 		u'  %s' %(s['label']),
# 		horizontalalignment='left',
# 		verticalalignment='center',
# 		color='k',
# 		alpha=.5
# 		)
# 
# 
# pylab.xlabel(u'δ47 (‰)')
# pylab.ylabel(u'Raw Δ47 (‰)')
# pylab.axis([xmin, xmax, None, None])
# 
# pylab.savefig('gas_correction.pdf')
