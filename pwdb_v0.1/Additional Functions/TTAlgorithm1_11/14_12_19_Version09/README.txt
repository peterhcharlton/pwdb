Copyright (c) 2012, Dr. Nicholas Gaddum
2012-Aug-10
Compiled in MATLAB 7.12.0 (R2011a)

Citations:

This work is the result of research documented in the following paper:

N R Gaddum, J Alastruey, P Beerbaum, P Chowienczyk, T Schaeffter,
A technical assessment of pulse wave velocity algorithms applied 
to non-invasive arterial waveforms, Annals of Biomedical Engineering,
2013, 41(12):2617-29. DOI: 10.1007/s10439-013-0854-y

Please cite this paper if this software is used for publication


Software Application:

Please read license for details of use.  Although use for pulse wave 
anlaysis this software should only be used for research application, 
and should not contribute to any clinical process including diagnosis,
assessment, patient planning or management.  No authors assume no
liabilities for the inaccuracies ofthe software.


Recomendations:

According to the paper noted above, the algorithm selection depends on
the type of waveform data, and the distance between the waveform 
measurements, (indicative of the similarity between wave shapes).  We 
recomend:
 - Low noise pressure data, (algorithm 1, Foot to foot (FTF))
 - Noisey, centrally measured data, (e.g. aortic Doppler/real time MRI),
	(algorithm 3, Least Squares)
 - Very noisey, centrally measured data, (e.g. aortic Doppler/real time 
	MRI), (algorithm 4, cross correlation of entire cycle)
 - Noisey data measured at a distance apart, (e.g. carotid/femoral artery
	tonometry), (algorithm 2, Foot to foot radius)


Instructions to RUN:

1. Read licence agreement.
2. It is assumed that the data consists of continuous waveform
	data, i.e. at least three cycles.  It is also assumed 
	that the waveform data was measured simultaneously, i.e.
	no temporal alignment is needed.
3. Install licenced MATLAB software on computer
4. Select TTALgorithm folder as the MATLAB current directory
5. Type:
	TT = TTAlgorithm(signal,f,algorithm,waveform,continuous,show)

	...to call function.

   where:

   signal = matrx consisting of two column vectors (side-by-side)
	with physiologcal waveform data at the two recorded sites
   f = sampling frequency (Hz)
   algorithm = 1, 2, 3, 4 to select foot to foot, foot to foot 
	radius, least squares or cross correlation of an entire
	waveform
   waveform = 1 or 2 to indicate whether the data is pressure/area
	or velocity/flow respectively
   continuous = 0, or 1, indicating a single wave (0) or multiple waves (1)
   show = output details/a figure showing the resulting located
	feet with which the calucation of TT data was facilitated

   TT = a column vector of transit times from the continous waveforms


Instructions to use Test Data:

High noise, low temporal resolution data (Doppler) and low noise, high 
temporal resolution data (Pressure) have been supplied to allow the 
user to test the solftware.

For the Doppler data use the following parameters:
signal		(included in supplied file)
f		(included in supplied file)
algorithm 	(1 to 4 depending on which algorithm to be applied)
waveform	(2, as this is velocity data)
continous	(1, as the data is continuous. 0, if a selection of 
			100 points (one cycle) is isused)
show		(1, to show the data with the located waveform features)
e.g.		TT = TTAlgorithm(signal,f,1,2,1,1)

For the Pressure data use the following parameters:
signal		(included in supplied file)
f		(included in supplied file)
algorithm 	(1 to 4 depending on which algorithm to be applied)
waveform	(1, as this is pressuredata)
continous	(1, as the data is continuous. 0, if a selection of 
			1000 points (one cycle) is isused)
show		(1, to show the data with the located waveform features)
e.g.		TT = TTAlgorithm(signal,f,1,1,1,1)


Trobleshooting:

1. Note that the wave 'foot' can be poorly lcoated if variable npoly_1
or npoly1_1 (in TTAlgorithm .m file) are poorly parameterised.  These
indicate the number of data points used within the sliding window to 
quantify the maximum systolic gradient, and the maximum end diastolic
gradient respectively.  Feel free to vary these two variables.


Revisions:
2012-Aug-06   Created function
2012-Aug-09   Removed replication of time parameters.  Only frequency 
	remaining.  Added pressure/area vs. velocity/flow parameter for
	Find_Landmarks observer.  Bug fixing.
2012-Aug-13   Revised README file.
2012-Aug-21   Pressure test data added to the Test Data File.  Added 
	recomendations.  Changes to the scanning window for algorithm 1.
2013-Aug-10   Bug fixing for high resolution pressure/flow data. Higher
	stability for foot location using the minimum radius
2014-Mar-24   Updated README file.  Added troubleshooting section
2014-Oct-02   Reduction of m files, and include option for single or
	multiple waveform analysis
2014-Dec-03   Bug fix for brachial-ankle, carotid-ankle data, where the
	feet have a greater spacing.  Cross correlation bug fix
2014-Dec-19   Removed and unused .mat file
2016-Jun-10   Bug fixing for low temporal resolution MRI flow data, and
	added signal observers for better stability of Cross Correlation