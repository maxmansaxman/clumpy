#! /usr/bin/env python
# -*- coding: utf-8 -*-
u'''
Web app to convert raw D47 values to "absolute" values
Mathieu Daëron (mathieu@daeron.fr), August 2014
'''

import xlrd
from flask import Flask, request, url_for, render_template
import numpy as np
from matplotlib import use
# use('Agg')
from matplotlib import rcParams
# rcParams['text.usetex'] = True
# rcParams['font.sans-serif'] = 'DejaVu Sans Mono'
import matplotlib.pyplot as mpl

# import matplotlib.font_manager
# print matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf')

def compute_CO2eqD47() :
	'''
	Theoretical equilibrium values for D47 in CO2 gas according to
	Wang et al. (2004) [http://dx.doi.org/10.1016/j.gca.2004.05.039],
	reported in the supplementary data of Dennis et al. (2011)
	[http://dx.doi.org/10.1016/j.gca.2011.09.025]
	'''

	eq_val = '''
	-83 1.8954\n-73 1.7530\n-63 1.6261\n-53 1.5126\n-43 1.4104\n-33 1.3182\n-23 1.2345
	-13 1.1584\n-3 1.0888\n7 1.0251\n17 0.9665\n27 0.9125\n37 0.8626\n47 0.8164\n57 0.7734
	67 0.7334\n87 0.6612\n97 0.6286\n107 0.5980\n117 0.5693\n127 0.5423\n137 0.5169
	147 0.4930\n157 0.4704\n167 0.4491\n177 0.4289\n187 0.4098\n197 0.3918\n207 0.3747
	217 0.3585\n227 0.3431\n237 0.3285\n247 0.3147\n257 0.3015\n267 0.2890\n277 0.2771
	287 0.2657\n297 0.2550\n307 0.2447\n317 0.2349\n327 0.2256\n337 0.2167\n347 0.2083
	357 0.2002\n367 0.1925\n377 0.1851\n387 0.1781\n397 0.1714\n407 0.1650\n417 0.1589
	427 0.1530\n437 0.1474\n447 0.1421\n457 0.1370\n467 0.1321\n477 0.1274\n487 0.1229
	497 0.1186\n507 0.1145\n517 0.1105\n527 0.1068\n537 0.1031\n547 0.0997\n557 0.0963
	567 0.0931\n577 0.0901\n587 0.0871\n597 0.0843\n607 0.0816\n617 0.0790\n627 0.0765
	637 0.0741\n647 0.0718\n657 0.0695\n667 0.0674\n677 0.0654\n687 0.0634\n697 0.0615
	707 0.0597\n717 0.0579\n727 0.0562\n737 0.0546\n747 0.0530\n757 0.0515\n767 0.0500
	777 0.0486\n787 0.0472\n797 0.0459\n807 0.0447\n817 0.0435\n827 0.0423\n837 0.0411
	847 0.0400\n857 0.0390\n867 0.0380\n877 0.0370\n887 0.0360\n897 0.0351\n907 0.0342
	917 0.0333\n927 0.0325\n937 0.0317\n947 0.0309\n957 0.0302\n967 0.0294\n977 0.0287
	987 0.0281\n997 0.0274\n1007 0.0268\n1017 0.0261\n1027 0.0255\n1037 0.0249\n1047 0.0244
	1057 0.0238\n1067 0.0233\n1077 0.0228\n1087 0.0223\n1097 0.0218
	'''

	eq_val = np.array([l.split() for l in eq_val.split('\n')[1:-1]], dtype='float')
	T = eq_val[:,0] - 0.15
	D47 = eq_val[:,1]

	# Quadratic polynomial fit as used by
	# Dennis et al. (2011) [http://dx.doi.org/10.1016/j.gca.2011.09.025]:
	eq_fit = np.polyfit( 1/(T+273.15), D47, 4 )
	return lambda T : \
		eq_fit[0] * (T+273.15)**-4 + \
		eq_fit[1] * (T+273.15)**-3 + \
		eq_fit[2] * (T+273.15)**-2 + \
		eq_fit[3] * (T+273.15)**-1 + \
		eq_fit[4]

CO2eqD47 = compute_CO2eqD47() # define the T law for equilibrium D47 values in CO2


def read_xls( content ) :
	'''
	reads the uploaded xls file and returns a list of measurements
	'''
	wb = xlrd.open_workbook( file_contents = content ) # open the contents of the uploaded xls file
	ws = wb.sheet_by_name( wb.sheet_names()[0] ) # get the first worksheet
	data = []
	for r in range( ws.nrows )[2:] : # for each row
		try :
			# read the cell values:
			label = str( ws.cell_value( r, 0) )
			d47 = float( ws.cell_value( r, 1) )
			rawD47 = float( ws.cell_value( r, 2) )
			srawD47 = float( ws.cell_value( r, 3) )
			data.append( {'label': label, 'd47': d47, 'rawD47': rawD47, 'srawD47': srawD47 } )
			if ws.cell_type( r, 4) == 2 :
				data[-1]['TCO2eq'] = float( ws.cell_value( r, 4) )
			if ws.cell_type( r, 5) == 2 :
				data[-1]['trueD47'] = float( ws.cell_value( r, 5) )
		except :
			pass
	for d in data :
		if 'TCO2eq' in d and 'trueD47' not in d :
			d['trueD47'] = CO2eqD47( d['TCO2eq'] ) # compute equilibrium value from gas T
	return data, hash( content )



def process_data( data, be_conservative = True ) :
	'''
	Computes the parameters for converting raw measurements
	into absolute reference frame values, using the following
	formulation:

		rawD47 = a * trueD47 + b * d47 + c

	The parameters (a,b,c) are computed using a LS fit of all
	the equilibrated gas and carbonate standard measurements.
	'''
	# design matrix:
	A = [ [ d['trueD47'] / d['srawD47'], d['d47'] / d['srawD47'], d['srawD47']**-1 ] for d in data if 'trueD47' in d ]
	# target values for the fit:
	Y = [ d['rawD47'] / d['srawD47' ] for d in data if 'trueD47' in d ]
	A,Y = np.array(A), np.array(Y)

	f = np.linalg.lstsq(A,Y.T)[0] # best-fit parameters
	CM = np.linalg.inv(np.dot(A.T,A)) # covariance matrix of fit parameters

	if be_conservative :
		# Scale up uncertainties in the fit parameters if the goodnes-of-fit is worse than average.
		# To some extent, this helps account for external errors in the gas line data.
		chi2 = sum( ( Y - np.dot( A, f ) )**2 )
		nf = len([ d for d in data if 'trueD47' in d ]) - 3
		if chi2 > nf :
			print "be_conservative: errors scaled by %.2f\n" %( chi2 / nf )
			CM = CM * chi2 / nf

	(a,b,c) = f
	for d in data :

		x = d['rawD47']
		y = d['d47']
		sx = d['srawD47']
		sy = 0.
		z = x / a - b/a * y - c/a
		d['corD47'] = z

		dzdx = a ** -1
		dzdy = -b / a
		dzda = -(x-b*y-c) * a ** -2
		dzdb = -y / a
		dzdc = -a ** -1

		v = np.array([ dzdx, dzdy, dzda, dzdb, dzdc ])
		C = np.zeros((5,5))
		C[0,0] = sx ** 2
		C[1,1] = sy ** 2
		C[2:,2:] = CM
		sz = (np.dot( v, np.dot( C, v.T ) ))**.5
		d['scorD47_all'] = sz

		C = np.zeros((5,5))
		C[0,0] = sx ** 2
		C[1,1] = sy ** 2
		sz = (np.dot( v, np.dot( C, v.T ) ))**.5
		d['scorD47_internal'] = sz

		C = np.zeros((5,5))
		C[2:,2:] = CM
		sz = (np.dot( v, np.dot( C, v.T ) ))**.5
		d['scorD47_model'] = sz

		# Although the above formulas hold true for linear combinations of normal variables
		# (Tellinghuisen, 2001) [http://dx.doi.org/10.1021/jp003484u],
		# Our Monte-Carlo simulations showed that distributions of D47 values
		# estimated using this code are undistinguishable from normal distributions,
		# so that 95% confidence limits correspond to +/- 2*sD47corrected.
		# These limits account for internal and external uncertainties on gas line data,
		# and for internal errors in sample measurements, but not for external uncertainties in the latter.
		# It is thus always recommended to check external reproducibility for
		# duplicate measurements of homogeneous samples.

	return f,CM

def plot_data( data, f, CM, filename ) :
	'''
	Generates a plot of the data overlaid with gas lines
	and model uncertainties.
	'''
	(a,b,c) = f
	fig = mpl.figure()
	fig.set_figwidth( 6 )
	fig.set_figheight( 6 )

	for d in data :
		msmarkers_kwargs = { 'ms':4, 'color':'w', 'mew':1 }
		errorbar_kwargs = { 'fmt':None, 'capsize':2, 'capthick':1, 'elinewidth':1 }
		if 'trueD47' in d :
			color = 'r'
			marker = 's'
			mpl.errorbar( d['d47'], d['rawD47'], 2*d['srawD47'], None, ecolor = color, **errorbar_kwargs )
			mpl.plot( d['d47'], d['rawD47'], marker, mec = color, **msmarkers_kwargs)
		else :
			color = 'b'
			marker = 'o'
			mpl.errorbar( d['d47'], d['rawD47'], 2*d['srawD47'], None, ecolor = color, **errorbar_kwargs )
			mpl.plot( d['d47'], d['rawD47'], marker, mec = color, **msmarkers_kwargs)

	xleft, xright, ybottom, ytop = mpl.axis()
	xi = np.linspace( xleft, xright )
	TCO2eq = list( set( [ d['TCO2eq'] for d in data if 'TCO2eq' in d ] ) )
	for T in TCO2eq :
		trueD47 = CO2eqD47( T )
		yi = a * trueD47 + b * xi + c
		dydx = b
		dyda = trueD47
		dydb = xi
		dydc = 1.
		C = CM
		syi = np.array( [ (np.dot( np.array([ dyda, xii, dydc ]), np.dot( C, np.array([ dyda, xii, dydc ]).T ) ))**.5 for xii in xi ] )
		mpl.fill_between( xi, yi-2*syi, yi+2*syi, color = 'r', alpha = 0.2 )

	mpl.axis( [ xleft, xright, None, None ] )
	xleft, xright, ybottom, ytop = mpl.axis()

	yi = np.arange( xleft, xright + 0.5, .5)
	xi = np.arange( ybottom, ytop + 0.05, .05)
	XI, YI = np.meshgrid( xi, yi)
	SI = XI * 0

	ny, nx = SI.shape
	for kx in range(nx) :
		for ky in range(ny) :
			x = XI[ky,kx]
			y = YI[ky,kx]

			dzda = -(x-b*y-c) * a ** -2
			dzdb = -y / a
			dzdc = -a ** -1

			v = np.array([ dzda, dzdb, dzdc ])
			C = CM
			SI[ky,kx] = (np.dot( v, np.dot( C, v.T ) ))**.5

	contour_kwargs = { 'ls':'-', 'c':'k', 'alpha':.25, 'zorder':-10 }
	# contour_label_kwargs = { 'horizontalalignment':'left', 'verticalalignment':'center', 'color':'cyan', 'alpha':1 }

	minS = SI.min()
	cs = mpl.contour( YI, XI, SI, [ np.ceil( minS * 1000 ) / 1000 + k * .002 for k in range(10) ], colors = contour_kwargs['c'], **contour_kwargs)
	mpl.clabel(cs)

	mpl.xlabel( ur'\u03B447 (\u2030)' )
	mpl.ylabel( ur'Raw \u039447 (\u2030)' )
	mpl.show()
	return
	# plot_path = '%d.jpg' % filename
	# if __name__ == '__main__' :
	# 	mpl.savefig( 'static/' + plot_path, dpi=150 )
	# else :
	# 	mpl.savefig( '/home/rambaldi/xlcor47/static/' + plot_path, dpi=150 )
	#
	# ### SAVE CSV ###
	# text_output = 'label	d47	rawD47	1SE	TCO2eq	absD47	corD47	SE(total)	SE(internal)	SE(model)\n\n'
	# for d in data :
	# 	if 'TCO2eq' in d :
	# 		text_output += '%s	%.4f	%.4f	%.4f	%.1f	%.4f	%.4f	%.4f	%.4f	%.4f\n' % ( d['label'], d['d47'], d['rawD47'], d['srawD47'], d['TCO2eq'], d['trueD47'], d['corD47'], d['scorD47_all'], d['scorD47_internal'], d['scorD47_model'] )
	# 	elif 'trueD47' in d :
	# 		text_output += '%s	%.4f	%.4f	%.4f		%.4f	%.4f	%.4f	%.4f	%.4f\n' % ( d['label'], d['d47'], d['rawD47'], d['srawD47'], d['trueD47'], d['corD47'], d['scorD47_all'], d['scorD47_internal'], d['scorD47_model'] )
	# 	else :
	# 		text_output += '%s	%.4f	%.4f	%.4f			%.4f	%.4f	%.4f	%.4f\n' % ( d['label'], d['d47'], d['rawD47'], d['srawD47'], d['corD47'], d['scorD47_all'], d['scorD47_internal'], d['scorD47_model'] )
	#
	# text_output += '\n\n'
	# text_output += 'MODEL PARAMETERS:\n'
	# text_output += '-----------------\n'
	# text_output += 'rawD47 = A * trueD47 + B * d47 + C\n'
	# text_output += 'A = %.6g\n' % a
	# text_output += 'B = %.6g\n' % b
	# text_output += 'C = %.6g\n' % c
	# text_output += 'Covariance matrix of best-fit (A,B,C) values:\n' % (-c/a)
	# text_output += '%s\n' % str( CM )
	# text_output += '---------------------------\n'
	# text_output += 'Scaling factor (A) = %.3g\n' % (a)
	# text_output += 'Heated gas slope (B/A) = %.3g\n' % (b/a)
	# text_output += 'Absolute D47 value of ref gas (-C/A) = %.3g permil\n' % (-c/a)
	#
	# return plot_path, text_output











app = Flask(__name__)

def allowed_file(filename):
	return '.' in filename and filename.rsplit('.', 1)[1] in ['xls']

@app.route('/', methods=['GET', 'POST'])
def upload_file():
	if request.method == 'POST':
		file = request.files['file']
		if file and allowed_file(file.filename):
			data, datahash = read_xls( file.read() )
			f,CM = process_data( data )
			plot_path, text_output = plot_data( data, f, CM, datahash )
			return render_template( 'output.html', plot_path = plot_path, lendata = '%d' %(len(data)+20), text_output = text_output )
	return render_template( 'input.html' )

if __name__ == '__main__' :
	app.debug = True
	app.run()
